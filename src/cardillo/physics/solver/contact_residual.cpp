#include "contact_residual.hpp"
#include "../../misc/contact_projection.hpp"
#include "../../misc/contact_rho.hpp"
#include "termination.hpp"

namespace cardillo::solver {

real_t contactNaturalMapResidual(physics::DynamicsAssembler& dyn, const VectorXr& u, real_t tolAbs, real_t tolRel) {
    const auto& W = dyn.W().asSparseRowMajor();
    const int Nc = (int)W.rows();
    if (Nc == 0) return 0;
    const int Nnfc = dyn.numFrictionlessContacts();
    const VectorXr& MinvDiag = dyn.MinvDiag();
    const VectorXr& bias = dyn.contactVVec();

    const VectorXr vPrev = W * dyn.vVec() + bias;
    const VectorXr xi = W * u + bias + dyn.restitutionVec().cwiseProduct(vPrev);

    VectorXr impulse = VectorXr::Zero(Nc);
    for (const auto& c : dyn.contacts()) {
        if (c.impulse_base_index < 0) continue;
        for (int k = 0; k < c.impulse_size; ++k) impulse[c.impulse_base_index + k] = c.last_impulse(k);
    }

    VectorXr rho = VectorXr::Zero(Nc);
    for (int row = 0; row < Nnfc; ++row) {
        real_t delassus = 0;
        for (Eigen::SparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(W, row); it; ++it) delassus += it.value() * it.value() * MinvDiag[it.col()];
        if (delassus > (real_t)0) rho[row] = (real_t)1 / delassus;
    }
    for (int row = Nnfc; row < Nc; row += 3) {
        const auto Wc = W.middleRows(row, 3);
        const Matrix33r G = Wc * MinvDiag.asDiagonal() * Wc.transpose();
        Vector3r contactRho = misc::computeContactRho(G, misc::ContactRhoStrategy::Split);
        if (!(contactRho.array() > (real_t)0).all()) contactRho = misc::computeContactRho(G, misc::ContactRhoStrategy::IsotropicDiagonal);
        rho.segment<3>(row) = contactRho;
    }

    const VectorXr& mu = dyn.muVec();
    VectorXr projected = impulse;
    for (int row = 0; row < Nnfc; ++row) projected[row] = misc::proxNormal(rho[row] * xi[row] - impulse[row]);
    for (int row = Nnfc; row < Nc; row += 3) {
        const Vector3r scaledVelocity = rho.segment<3>(row).cwiseProduct(xi.segment<3>(row));
        misc::proxCoulomb(scaledVelocity - impulse.segment<3>(row), mu[row], projected.segment<3>(row));
    }

    return termination::naturalMapError(impulse, projected, rho, xi, xi, Nnfc, tolAbs, tolRel);
}

}  // namespace cardillo::solver
