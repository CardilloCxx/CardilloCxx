#include "local_contact_newton.hpp"

#include <Eigen/LU>
#include <cmath>

namespace cardillo::solver {

namespace {

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
    if (G(0, 0) <= eps) return false;
    const real_t rhoN = (real_t)1 / G(0, 0);

    const real_t a = G(1, 1), d = G(2, 2), b = G(1, 2);
    const real_t disc = std::max((real_t)0, (a - d) * (a - d) + (real_t)4 * b * b);
    const real_t lambdaMaxTT = (real_t)0.5 * ((a + d) + std::sqrt(disc));
    if (lambdaMaxTT <= eps) return false;
    const real_t rhoT = (real_t)1 / lambdaMaxTT;

    const Vector3r rho(rhoN, rhoT, rhoT);
    Matrix33r M = Matrix33r::Identity() - rho.asDiagonal() * G;
    const Vector3r q = r + G * lambda;
    const Vector3r c = rho.asDiagonal() * q;

    auto phiAt = [&](const Vector3r& lam) {
        const Vector3r y = M * lam + c;
        const CaseEval ce = evalCase(M, y, mu);
        return std::make_pair(Vector3r(lam - ce.g), Matrix33r(Matrix33r::Identity() - ce.Jg));
    };

    Vector3r lamIter = lambda;
    auto [phi, Jphi] = phiAt(lamIter);

    for (int it = 0; it < params.maxIters; ++it) {
        if (phi.norm() < params.tol) {
            lambda = lamIter;
            return true;
        }

        Eigen::FullPivLU<Matrix33r> lu(Jphi);
        if (!lu.isInvertible()) return false;
        const Vector3r delta = lu.solve(-phi);
        if (!delta.allFinite()) return false;

        const Vector3r lamTrial = lamIter + delta;
        auto [phiTrial, JphiTrial] = phiAt(lamTrial);
        if (it > 0 && phiTrial.norm() > phi.norm()) return false;

        lamIter = lamTrial;
        phi = phiTrial;
        Jphi = JphiTrial;
    }

    return false;  // maxIters exhausted without meeting tol
}

}  // namespace cardillo::solver
