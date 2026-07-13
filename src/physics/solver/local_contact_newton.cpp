#include "local_contact_newton.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

namespace cardillo::solver {

namespace {

// NewtonRhoStrategy::FullSpectral (see local_contact_newton.hpp): rho = 1/lambda_max(sym(G)),
// applied to all three components -- matches Siconos's compute_rho_spectral_norm() exactly
// (fc3d_AlartCurnier_functions.c, via op3x3.h's eig_3x3(), documented there as "eigenvalues of a
// SYMMETRIC 3x3 real matrix"). G need not be exactly symmetric in this codebase once an
// implicit-gyroscopic body is involved (see condensed_assembler.cpp) -- symmetrizing first, rather
// than using a general (possibly-complex-eigenvalue) solver, mirrors what Siconos's own
// closed-form routine assumes and keeps this a real scalar step size in every case.
real_t rhoFullSpectral(const Matrix33r& G, real_t eps) {
    const Matrix33r Gsym = (real_t)0.5 * (G + G.transpose());
    Eigen::SelfAdjointEigenSolver<Matrix33r> es(Gsym, Eigen::EigenvaluesOnly);
    const real_t lambdaMax = es.eigenvalues()(2);  // ascending order -- largest is last
    return (lambdaMax > eps) ? (real_t)1 / lambdaMax : (real_t)-1;
}

// One case-split evaluation of the Alart-Curnier normal map, in PGS's lambda[0]<=0 convention.
// Given the (per-Newton-call-constant) affine trial map y(lambda) = M*lambda + c, computes the
// projected target g(lambda), its Jacobian dg/dlambda, and Phi = lambda - g.
//
// Derivation note: y(lambda') is built from the LOCAL linear model of how this block's own
// residual r(lambda') changes as ONLY this block's impulse varies (all other blocks' impulses,
// hence u_corr's other contributions, held fixed): r(lambda') = r_current - G*(lambda'-lambda_cur)
// (the minus sign because increasing this block's own lambda increases u_corr via
// `u_corr += Minv*J^T*Delta_lambda`, which in turn DECREASES r = rhs - matvec(u_corr) -- this is
// this codebase's own u_corr/rhs convention, not a generic textbook one, and was checked directly
// against CondensedSolver's own linear/projection update `lambda_new = lambda_cur + GiiInv*r`,
// which is exactly the lambda'=M^{-1}... fixed point of r(lambda')=0 under this sign). The
// resulting one-step trial state (before projection) is y(lambda') = lambda' + rho*r(lambda'),
// with q := r_current + G*lambda_cur, M := I - rho*G, c := rho*q -- so that y(lambda_cur) reduces
// to exactly lambda_cur + rho*r_current, matching the plain relaxation*GiiInv*r update's role as
// "trial state before projection" for rho == GiiInv.
struct CaseEval {
    Vector3r g;
    Matrix33r Jg;  // dg/dlambda
};

CaseEval evalCase(const Matrix33r& M, const Vector3r& y, real_t mu) {
    CaseEval out;
    const real_t yN = y[0];
    const Eigen::Matrix<real_t, 2, 1> yT = y.template segment<2>(1);
    const real_t s = yT.norm();

    if (yN >= (real_t)0) {
        // Separation: no contact force.
        out.g.setZero();
        out.Jg.setZero();
        return out;
    }

    if (s <= -mu * yN) {
        // Sticking: unconstrained trial state is already inside (or on) the friction cone.
        out.g = y;
        out.Jg = M;
        return out;
    }

    // Sliding: project the tangential trial state onto the disk of radius -mu*yN.
    const Eigen::Matrix<real_t, 2, 1> that_hat = yT / s;
    out.g[0] = yN;
    out.g.template segment<2>(1) = -mu * yN * that_hat;

    const Eigen::Matrix<real_t, 1, 3> aN = M.row(0);
    const Eigen::Matrix<real_t, 2, 3> BT = M.block<2, 3>(1, 0);
    const Eigen::Matrix<real_t, 2, 2> proj = Eigen::Matrix<real_t, 2, 2>::Identity() - that_hat * that_hat.transpose();
    const Eigen::Matrix<real_t, 2, 3> dgT = -mu * that_hat * aN - (mu * yN / s) * proj * BT;

    out.Jg.row(0) = aN;
    out.Jg.block<2, 3>(1, 0) = dgT;
    return out;
}

}  // namespace

bool solveContactBlockNewtonAC(const Matrix33r& G, const Vector3r& r, real_t mu, Vector3r& lambda, const NewtonACParams& params) {
    const real_t eps = (real_t)1e-12;

    Vector3r rho;
    if (params.rhoStrategy == NewtonRhoStrategy::FullSpectral) {
        // see Acary2018 equation (110): https://doi.org/10.1007/978-3-319-75972-2_10
        const real_t rhoAll = rhoFullSpectral(G, eps);
        if (rhoAll <= (real_t)0) return false;
        rho = Vector3r(rhoAll, rhoAll, rhoAll);
    } else {
        // see Acary2018 equation (111): https://doi.org/10.1007/978-3-319-75972-2_10
        if (G(0, 0) <= eps) return false;
        const real_t rhoN = (real_t)1 / G(0, 0);

        // Closed-form (quadratic-formula) largest eigenvalue of the symmetrized 2x2 tangential
        // block -- already exact for a symmetric input, so routing this through a generic
        // SelfAdjointEigenSolver would add setup overhead for no accuracy/robustness gain. `b` is
        // averaged from both off-diagonal entries (not just G(1,2)) because G need not be
        // symmetric here: an implicit-gyroscopic body (moreau.implicit_gyroscopy=true, see
        // condensed_assembler.cpp's gyroMinvBlocks) makes Gii -- and hence this G -- genuinely
        // non-symmetric, a live, reachable configuration, not a hypothetical one. Using only
        // G(1,2) would silently ignore G(2,1) instead of accounting for it. Symmetrizing first
        // also matches rhoFullSpectral's convention and, crucially, guarantees the discriminant
        // below is exactly the real quantity (a-d)^2+4b^2 >= 0 for the symmetrized block -- the
        // eigenvalues of the ACTUAL (unsymmetrized, possibly non-normal) 2x2 block could otherwise
        // be genuinely complex (e.g. a rotation-like block), which this real closed-form formula
        // cannot represent.
        const real_t a = G(1, 1), d = G(2, 2), b = (real_t)0.5 * (G(1, 2) + G(2, 1));
        const real_t disc = std::max((real_t)0, (a - d) * (a - d) + (real_t)4 * b * b);
        const real_t lambdaMaxTT = (real_t)0.5 * ((a + d) + std::sqrt(disc));
        if (lambdaMaxTT <= eps) return false;
        const real_t rhoT = (real_t)1 / lambdaMaxTT;

        rho = Vector3r(rhoN, rhoT, rhoT);
    }
    Matrix33r M = Matrix33r::Identity() - rho.asDiagonal() * G;
    const Vector3r q = r + G * lambda;
    const Vector3r c = rho.asDiagonal() * q;

    // Writes into out-params rather than returning a fresh pair -- every call site below reuses
    // the same pre-declared storage instead of constructing new Vector3r/Matrix33r each time.
    auto phiAt = [&](const Vector3r& lam, Vector3r& phiOut, Matrix33r& JphiOut) {
        const Vector3r y = M * lam + c;
        const CaseEval ce = evalCase(M, y, mu);
        phiOut = lam - ce.g;
        JphiOut = Matrix33r::Identity() - ce.Jg;
    };

    Vector3r lamIter = lambda;
    Vector3r phi, delta, lamTrial, phiTrial;
    Matrix33r Jphi, JphiInv, JphiTrial;
    bool invertible;
    phiAt(lamIter, phi, Jphi);

    for (int it = 0; it < params.maxIters; ++it) {
        if (phi.norm() < params.tol) {
            lambda = lamIter;
            return true;
        }

        // computeInverseWithCheck's invertibility flag falls out of computing the inverse itself
        // (no extra decomposition/cost) and catches the exact/near-exact singular case cleanly;
        // delta.allFinite() below is a second, independent guard against a bad step slipping
        // through (see CONDENSED_SOLVER_REPORT.md for the measured gap and why it's latent here).
        Jphi.computeInverseWithCheck(JphiInv, invertible);
        if (!invertible) return false;
        delta = JphiInv * (-phi);
        if (!delta.allFinite()) return false;

        lamTrial = lamIter + delta;
        phiAt(lamTrial, phiTrial, JphiTrial);
        if (it > 0 && phiTrial.norm() > phi.norm()) return false;

        lamIter = lamTrial;
        phi = phiTrial;
        Jphi = JphiTrial;
    }

    return false;  // maxIters exhausted without meeting tol
}

}  // namespace cardillo::solver
