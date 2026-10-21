#include "local_contact_newton.hpp"

#include <cmath>

namespace cardillo::solver {

namespace {

// One case-split evaluation of the Alart-Curnier map in the negative variant of the theory manual
// (eq:newton_trial, eq:alart_curnier, same as misc/contact_projection.hpp): the trial point
//   y(lambda') = rho*xi(lambda') - lambda' = -M*lambda' + c,   M = I - rho*G,   c = rho*(xi_cur - G*lambda_cur),
// follows from the local affine model xi(lambda') = xi_cur + G*(lambda' - lambda_cur) of this block's
// row velocity with all other blocks frozen, and p(y) = -proj_C(y) is the new percussion. Computes
// p and dp/dlambda = dp/dy * (-M); Phi = lambda - p.
struct CaseEval {
    Vector3r g;    // p(y)
    Matrix33r Jg;  // dp/dlambda
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
        // Sticking: -y is already inside (or on) the friction disk, p = -y, dp/dlambda = M.
        out.g = -y;
        out.Jg = M;
        return out;
    }

    // Sliding: p = (-yN, -proj_D(yT)) with the disk radius -mu*yN, i.e. p_T = mu*yN*t_hat.
    const Eigen::Matrix<real_t, 2, 1> that_hat = yT / s;
    out.g[0] = -yN;
    out.g.template segment<2>(1) = mu * yN * that_hat;

    const Eigen::Matrix<real_t, 1, 3> aN = M.row(0);
    const Eigen::Matrix<real_t, 2, 3> BT = M.block<2, 3>(1, 0);
    const Eigen::Matrix<real_t, 2, 2> proj = Eigen::Matrix<real_t, 2, 2>::Identity() - that_hat * that_hat.transpose();
    const Eigen::Matrix<real_t, 2, 3> dgT = -mu * that_hat * aN - (mu * yN / s) * proj * BT;

    out.Jg.row(0) = aN;
    out.Jg.block<2, 3>(1, 0) = dgT;
    return out;
}

}  // namespace

bool solveContactBlockNewtonAC(const Matrix33r& G, const Vector3r& xi, real_t mu, Vector3r& lambda, const NewtonACParams& params) {
    const real_t eps = (real_t)1e-12;

    const Vector3r rho = misc::computeContactRho(G, params.rhoStrategy, eps);
    if ((rho.array() <= (real_t)0).any()) return false;

    Matrix33r M = Matrix33r::Identity() - rho.asDiagonal() * G;
    const Vector3r q = xi - G * lambda;
    const Vector3r c = rho.asDiagonal() * q;

    // Writes into out-params rather than returning a fresh pair -- every call site below reuses
    // the same pre-declared storage instead of constructing new Vector3r/Matrix33r each time.
    auto phiAt = [&](const Vector3r& lam, Vector3r& phiOut, Matrix33r& JphiOut) {
        const Vector3r y = c - M * lam;
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
        // through.
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
