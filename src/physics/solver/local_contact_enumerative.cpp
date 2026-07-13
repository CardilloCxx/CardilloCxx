#include "local_contact_enumerative.hpp"

#include "quartic_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cardillo::solver {

namespace {

// One accepted candidate before the (rare, non-symmetric-G-only) tie-break -- see this file's
// header comment on non-uniqueness.
struct Candidate {
    Vector3r lambda;
    bool valid{false};
};

// Case-agnostic Alart-Curnier-style fixed-point residual (rho=1), used ONLY to rank multiple
// simultaneously-valid candidates against each other. Phi(lambda) = lambda - proj_disk(lambda +
// u(lambda)) is exactly zero iff `lambda` solves the local NCP, for ANY rho > 0 -- a fixed rho=1 is
// enough for this coarse ranking purpose; there's no need to duplicate local_contact_newton.cpp's
// own per-block-optimized rho machinery just to break a tie here.
real_t candidateResidual(const Matrix33r& G, const Vector3r& q, real_t mu, const Vector3r& lam) {
    const Vector3r u = q - G * lam;
    const Vector3r y = lam + u;
    const real_t yN = y[0];
    const Eigen::Matrix<real_t, 2, 1> yT = y.template segment<2>(1);

    Vector3r g = Vector3r::Zero();
    if (yN < (real_t)0) {
        const real_t s = yT.norm();
        if (s <= -mu * yN) {
            g = y;
        } else {
            g[0] = yN;
            g.template segment<2>(1) = (-mu * yN / s) * yT;
        }
    }
    return (lam - g).norm();
}

}  // namespace

bool solveContactBlockEnumerative(const Matrix33r& G, const Vector3r& r, real_t mu, Vector3r& lambda, const EnumerativeParams& params) {
    const real_t eps = params.eps;

    // `r` is a residual relative to the INCOMING lambda (matches local_contact_newton.cpp:118's own
    // q = r + G*lambda exactly) -- this solver produces an ABSOLUTE new lambda, so it needs the
    // true free term q, not raw r. See this file's header comment.
    const Vector3r q = r + G * lambda;

    // Verified sign/notation translation from the paper's r_N_paper>=0 convention to this
    // codebase's lambda[0]<=0 convention (see CONDENSED_SOLVER_REPORT.md and the planning record
    // for the numeric verification against planted analytic solutions, both symmetric and
    // non-symmetric -- an earlier hand-derivation of these signs was WRONG and is not what's here).
    const real_t wN = G(0, 0);
    const Eigen::Matrix<real_t, 1, 2> wTr = -G.template block<1, 2>(0, 1);
    const Eigen::Matrix<real_t, 2, 1> wTc = G.template block<2, 1>(1, 0);
    const Eigen::Matrix<real_t, 2, 2> WT = -G.template block<2, 2>(1, 1);
    const real_t bN = q(0);
    const Eigen::Matrix<real_t, 2, 1> bT = q.template segment<2>(1);

    Candidate best;
    real_t bestResidual = std::numeric_limits<real_t>::infinity();
    auto consider = [&](const Candidate& c) {
        if (!c.valid) return;
        const real_t res = candidateResidual(G, q, mu, c.lambda);
        if (res < bestResidual) {
            bestResidual = res;
            best = c;
        }
    };

    // Take-off: lambda=0, valid iff the free-term normal velocity is already non-negative.
    consider(Candidate{Vector3r::Zero(), bN >= -eps});

    // Sticking: the unique lambda with u(lambda)=0, i.e. lambda = G^{-1}*q -- a fresh 3x3 inverse
    // (matching how local_contact_newton.cpp inverts its own Jacobian), not blk.GiiInv, which isn't
    // reachable from this function's signature. Valid iff inside (or on) the friction cone.
    {
        Matrix33r Ginv;
        bool invertible;
        G.computeInverseWithCheck(Ginv, invertible);
        if (invertible) {
            const Vector3r lamStick = Ginv * q;
            const real_t s = lamStick.template segment<2>(1).norm();
            const bool ok = lamStick[0] <= eps && s <= -mu * lamStick[0] + eps;
            Vector3r lamClamped = lamStick;
            lamClamped[0] = std::min(lamClamped[0], (real_t)0);
            consider(Candidate{lamClamped, ok});
        }
    }

    if (mu <= eps) {
        // Frictionless: the sliding case collapses to the plain normal LCP, r_T=0, solved directly
        // (Bonnefon & Daviet's own Appendix B.2 special-cases this too -- the general sliding
        // machinery below is genuinely uninformative at mu=0, not just slower: with r_T fixed at 0,
        // the cone-boundary condition ||r_T||=mu*r_N is trivially satisfied for ANY r_N, so there is
        // no scalar equation left to solve for alpha).
        Vector3r lamFrictionless = Vector3r::Zero();
        if (std::abs(wN) > eps) {
            lamFrictionless[0] = std::min(bN / wN, (real_t)0);
            consider(Candidate{lamFrictionless, true});
        }
    } else if (std::abs(wN) > eps) {
        // Sliding: reduces to a scalar quartic in alpha (see local_contact_enumerative.hpp and
        // CONDENSED_SOLVER_REPORT.md for the derivation, verified both symbolically and against
        // planted analytic solutions for non-symmetric G).
        const Eigen::Matrix<real_t, 2, 2> Wbar = WT - (wTc * wTr) / wN;
        const Eigen::Matrix<real_t, 2, 1> bbar = bT - (bN / wN) * wTc;

        const real_t a1 = Wbar(0, 0) + Wbar(1, 1);
        const real_t a0 = Wbar(0, 0) * Wbar(1, 1) - Wbar(0, 1) * Wbar(1, 0);
        Eigen::Matrix<real_t, 2, 2> adjWbar;
        adjWbar << Wbar(1, 1), -Wbar(0, 1), -Wbar(1, 0), Wbar(0, 0);
        const Eigen::Matrix<real_t, 2, 1> B = adjWbar * bbar;

        const real_t n2 = bbar.dot(bbar), n1 = (real_t)2 * B.dot(bbar), n0 = B.dot(B);
        const real_t q2 = -bN, q1 = wTr.dot(bbar) - bN * a1, q0 = wTr.dot(B) - bN * a0;

        const real_t a4 = -mu * mu * q2 * q2;
        const real_t a3 = (real_t)-2 * mu * mu * q2 * q1;
        const real_t a2q = wN * wN * n2 - mu * mu * ((real_t)2 * q2 * q0 + q1 * q1);
        const real_t a1q = wN * wN * n1 - (real_t)2 * mu * mu * q1 * q0;
        const real_t a0q = wN * wN * n0 - mu * mu * q0 * q0;

        real_t roots[4];
        const int n = solveQuarticReal(a4, a3, a2q, a1q, a0q, roots);
        for (int i = 0; i < n; ++i) {
            const real_t alpha = roots[i];
            if (alpha <= eps) continue;
            const real_t Delta = alpha * alpha + a1 * alpha + a0;
            if (std::abs(Delta) <= eps) continue;

            const Eigen::Matrix<real_t, 2, 1> rTpaper = -(B + alpha * bbar) / Delta;
            const real_t rNpaper = (q2 * alpha * alpha + q1 * alpha + q0) / (wN * Delta);
            if (rNpaper < -eps) continue;  // squaring the cone-magnitude equation lost this sign -- check explicitly

            Vector3r lamSlide;
            lamSlide[0] = std::min(-rNpaper, (real_t)0);
            lamSlide.template segment<2>(1) = rTpaper;
            consider(Candidate{lamSlide, true});
        }
    }

    if (!best.valid) return false;
    lambda = best.lambda;
    return true;
}

}  // namespace cardillo::solver
