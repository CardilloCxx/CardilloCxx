#include "natural_map.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <Eigen/Eigenvalues>

#include "../../misc/contact_projection.hpp"
#include "../assembly/dynamics_assembler.hpp"

namespace cardillo::solver {

namespace {

// Entry (i, j) of W M^-1 W^T with the diagonal mass, from the sorted rows i and j of W.
real_t delassusEntry(const CsrMatrix& W, int i, int j, const VectorXr& MinvDiag) {
    real_t sum = 0;
    CsrMatrix::InnerIterator a(W, i), b(W, j);
    while (a && b) {
        if (a.col() < b.col()) {
            ++a;
        } else if (b.col() < a.col()) {
            ++b;
        } else {
            sum += a.value() * b.value() * MinvDiag[a.col()];
            ++a;
            ++b;
        }
    }
    return sum;
}

}  // namespace

ContactNaturalMap::ContactNaturalMap(const CsrMatrix& W, const VectorXr& MinvDiag, const VectorXr& mu, int numFrictionlessRows)
    : m_numFrictionlessRows(numFrictionlessRows), m_mu(mu), m_rho(VectorXr::Zero(W.rows())), m_invRho(VectorXr::Zero(W.rows())) {
    for (int row = 0; row < numFrictionlessRows; ++row) {
        const real_t G = delassusEntry(W, row, row, MinvDiag);
        if (G > (real_t)0) {
            m_rho[row] = (real_t)1 / G;
            m_invRho[row] = G;
        }
    }
    for (int row = numFrictionlessRows; row + 3 <= (int)W.rows(); row += 3) {
        Matrix33r G;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) G(i, j) = G(j, i) = delassusEntry(W, row + i, row + j, MinvDiag);
        Eigen::SelfAdjointEigenSolver<Matrix33r> eigen;
        eigen.computeDirect(G, Eigen::EigenvaluesOnly);
        const real_t lambdaMax = eigen.eigenvalues()[2];
        if (lambdaMax > (real_t)0) {
            m_rho.segment<3>(row).setConstant((real_t)1 / lambdaMax);
            m_invRho.segment<3>(row).setConstant(lambdaMax);
        }
    }
}

real_t ContactNaturalMap::residual(const VectorXr& lambda, const VectorXr& xi) const {
    return maxResidual(lambda, xi, nullptr);
}

real_t ContactNaturalMap::residual(const VectorXr& lambda, const VectorXr& xi, const VectorXr& frozenShift) const {
    return maxResidual(lambda, xi, &frozenShift);
}

real_t ContactNaturalMap::maxResidual(const VectorXr& lambda, const VectorXr& xi, const VectorXr* frozenShift) const {
    // Squares, with one square root at the end.
    const int Nnfc = m_numFrictionlessRows;
    const int numFrictional = ((int)m_rho.size() - Nnfc) / 3;
    real_t residual2 = 0;
    bool finite = true;

#pragma omp parallel for schedule(static) reduction(max : residual2) reduction(&& : finite)
    for (int row = 0; row < Nnfc; ++row) {
        const real_t rho = m_rho[row];
        if (!(rho > (real_t)0)) continue;
        const real_t r = (lambda[row] - std::max((real_t)0, lambda[row] - rho * xi[row])) * m_invRho[row];
        finite = finite && std::isfinite(r);
        residual2 = std::max(residual2, r * r);
    }

#pragma omp parallel for schedule(static) reduction(max : residual2) reduction(&& : finite)
    for (int c = 0; c < numFrictional; ++c) {
        const int row = Nnfc + 3 * c;
        const real_t rho = m_rho[row];
        if (!(rho > (real_t)0)) continue;
        const Vector3r lam = lambda.segment<3>(row);
        Vector3r corrected = xi.segment<3>(row);
        corrected[0] += frozenShift ? (*frozenShift)[row] : misc::deSaxce(corrected, m_mu[row]);
        const Vector3r r = (lam - misc::projectCoulombCone(lam - rho * corrected, m_mu[row])) * m_invRho[row];
        const real_t rN2 = r[0] * r[0];
        const real_t rT2 = r[1] * r[1] + r[2] * r[2];
        finite = finite && std::isfinite(rN2 + rT2);
        residual2 = std::max(residual2, std::max(rN2, rT2));
    }

    return finite ? std::sqrt(residual2) : std::numeric_limits<real_t>::quiet_NaN();
}

StepResidual stepResidual(physics::DynamicsAssembler& dyn, const VectorXr& u, ContactLaw law, real_t dt, real_t theta) {
    StepResidual result;

    const int Nc = dyn.numContactRows();
    if (Nc > 0) {
        const CsrMatrix& W = dyn.W().asSparseRowMajor();
        const int Nnfc = dyn.numFrictionlessContacts();
        const VectorXr gammaPrev = W * dyn.vVec() + dyn.contactVVec();  // pre-step contact velocity
        const VectorXr xi = W * u + dyn.contactVVec() + dyn.restitutionVec().cwiseProduct(gammaPrev);
        const VectorXr lambda = dyn.contactPercussions();

        const ContactNaturalMap map(W, dyn.MinvDiag(), dyn.muVec(), Nnfc);
        // Coulomb's law, the physics error, for every solver. The solver's own law is the same, except
        // for the convexified law of the conic backends, which freezes the shift at mu ||gamma_T^-||
        // (conic/conic_scaling.hpp); the difference of the two is the error of the convexification.
        result.coulomb = map.residual(lambda, xi);
        if (law == ContactLaw::Coulomb) {
            result.contact = result.coulomb;
        } else {
            VectorXr shift = VectorXr::Zero(Nc);
            for (int row = Nnfc; row + 3 <= Nc; row += 3) shift[row] = misc::deSaxce(gammaPrev.segment<3>(row), dyn.muVec()[row]);
            result.contact = map.residual(lambda, xi, shift);
        }
    }

    // Spring and damper rows: W_i^T u + b_i + c_i Lambda_i = 0 (DynamicsAssembler::springBias()).
    auto rowResidual = [&](const VectorXr& rows) {
        if (rows.size() > 0) result.bilateral = std::max(result.bilateral, rows.cwiseAbs().maxCoeff());
    };
    const int Ns = dyn.numSprings();
    if (Ns > 0 && dyn.Lambda_g().size() == Ns)
        rowResidual(dyn.Wg() * u + dyn.springBias(dt, theta) + (physics::DynamicsAssembler::springComplianceScale(dt, theta) * dyn.Cdiag()).cwiseProduct(dyn.Lambda_g()));
    const int Nd = dyn.numDampers();
    if (Nd > 0 && dyn.Lambda_gamma().size() == Nd)
        rowResidual(dyn.Wgamma() * u + dyn.damperBias(dt, theta) + (physics::DynamicsAssembler::damperComplianceScale(dt, theta) * dyn.Adiag()).cwiseProduct(dyn.Lambda_gamma()));

    return result;
}

}  // namespace cardillo::solver
