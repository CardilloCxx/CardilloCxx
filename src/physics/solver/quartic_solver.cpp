#include "quartic_solver.hpp"

#include <unsupported/Eigen/Polynomials>

#include <algorithm>
#include <cmath>
#include <limits>

namespace cardillo::solver {

namespace {

// Combined Horner scheme (synthetic division): one pass evaluates both P(x) and P'(x).
void hornerEval(real_t a4, real_t a3, real_t a2, real_t a1, real_t a0, real_t x, real_t& P, real_t& dP) {
    P = a4;
    dP = (real_t)0;
    for (real_t c : {a3, a2, a1, a0}) {
        dP = dP * x + P;
        P = P * x + c;
    }
}

// Up to 2 Newton corrections against the ORIGINAL (un-normalized) quartic. Companion-matrix
// eigenvalues are already close to machine precision for well-scaled problems, so this is a cheap
// accuracy touch-up, not a from-scratch iteration -- guarded against a non-finite/zero derivative
// and against a step that makes the residual worse (mirrors local_contact_newton.cpp's own
// non-increase guard style).
real_t polishRoot(real_t a4, real_t a3, real_t a2, real_t a1, real_t a0, real_t x0) {
    real_t x = x0, P, dP;
    hornerEval(a4, a3, a2, a1, a0, x, P, dP);
    for (int it = 0; it < 2; ++it) {
        if (!std::isfinite(dP) || dP == (real_t)0) break;
        const real_t xTrial = x - P / dP;
        if (!std::isfinite(xTrial)) break;
        real_t Ptrial, dPtrial;
        hornerEval(a4, a3, a2, a1, a0, xTrial, Ptrial, dPtrial);
        if (std::abs(Ptrial) > std::abs(P)) break;
        x = xTrial;
        P = Ptrial;
        dP = dPtrial;
    }
    return x;
}

}  // namespace

int solveQuarticReal(real_t a4, real_t a3, real_t a2, real_t a1, real_t a0, real_t roots[4]) {
    // Deliberately close to machine epsilon, not a generous "looks negligible" threshold: the
    // companion-matrix solve below (Eigen::PolynomialSolver) handles a small-but-present leading
    // coefficient correctly on its own (verified down to a4/scale ~ 1e-11 against numpy.roots, see
    // CONDENSED_SOLVER_REPORT.md) -- a looser threshold here would silently DROP a genuine
    // (possibly huge-magnitude) root by mistaking "small" for "zero", which is worse than passing a
    // slightly-imprecise near-zero coefficient down to a solver that already copes with it.
    const real_t scale = std::max({std::abs(a4), std::abs(a3), std::abs(a2), std::abs(a1), std::abs(a0), (real_t)1});
    const real_t eps = (real_t)4 * std::numeric_limits<real_t>::epsilon() * scale;

    // Ascending order, as Eigen::PolynomialSolver expects. Walk the degree down while the current
    // leading coefficient is negligible -- see quartic_solver.hpp's comment on why a4~=0 is a real,
    // expected case for this solver's caller, not a corner case to guard against and forget.
    const real_t coeffsAsc[5] = {a0, a1, a2, a3, a4};
    int degree = 4;
    while (degree > 0 && std::abs(coeffsAsc[degree]) <= eps) --degree;

    int count = 0;
    if (degree == 0) {
        return 0;  // either "0=0" (no information) or a nonzero constant (no roots) -- both give 0
    }

    VectorXr c(degree + 1);
    for (int i = 0; i <= degree; ++i) c(i) = coeffsAsc[i];
    Eigen::PolynomialSolver<real_t, Eigen::Dynamic> solver;
    solver.compute(c);

    const auto& polyRoots = solver.roots();
    for (Eigen::Index i = 0; i < polyRoots.size(); ++i) {
        const real_t re = polyRoots[i].real();
        const real_t im = polyRoots[i].imag();
        // A genuine real (possibly repeated) root of a real polynomial can still come out of the
        // companion-matrix eigendecomposition with a small nonzero imaginary part (repeated roots
        // are the classic hard case for eigenvalue-based root finders) -- classify by size relative
        // to the root's own magnitude, not the polynomial's coefficient scale.
        const real_t rootScale = std::max(std::abs(re), (real_t)1);
        if (std::abs(im) <= (real_t)1e-6 * rootScale) {
            roots[count++] = re;
        }
    }

    for (int i = 0; i < count; ++i) roots[i] = polishRoot(a4, a3, a2, a1, a0, roots[i]);
    return count;
}

}  // namespace cardillo::solver
