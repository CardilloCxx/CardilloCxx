#pragma once

#include "../../misc/types.hpp"

namespace cardillo::solver {

// Robust real-root finder for a general quartic a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0, given as
// five plain coefficients (no polynomial/matrix wrapper type -- this is a standalone numerical
// utility with zero contact-specific knowledge, reusable anywhere a quartic needs real roots).
//
// Implementation: Eigen's unsupported companion-matrix polynomial solver (a balanced real Schur
// eigendecomposition of the (degree x degree) companion matrix), each real root then polished with
// up to 2 Newton iterations (Horner evaluation) against the ORIGINAL, un-normalized coefficients.
//
// A hand-rolled closed-form (Ferrari + trigonometric resolvent-cubic) implementation was tried
// first for speed, and validated against thousands of random + adversarial quartics (repeated
// roots, a near-zero leading coefficient at several magnitudes, widely-separated root magnitudes --
// exactly the categories this solver's actual caller can hit: local_contact_enumerative.cpp's
// sliding case has a4 = -mu^2*bN^2, which genuinely vanishes as bN -> 0, a real physical case, not
// a synthetic corner case). That closed-form attempt FAILED two of those categories outright -- it
// silently dropped genuine double roots entirely, and returned roots wrong by 10+ orders of
// magnitude whenever the leading coefficient was small-but-nonzero relative to the rest (b3=a3/a4
// blows up under monic normalization, and the depression formulas' b3^4 term then swamps the
// O(1) terms it needs to retain). Eigen::PolynomialSolver passed every one of those same cases
// exactly (cross-checked against numpy.roots as an independent reference). See
// CONDENSED_SOLVER_REPORT.md for the numbers. A 4x4 (or smaller, see the degree-detection below)
// real eigenvalue problem is still cheap -- this is not a large iterative solve.
//
// Degree-detection: if |a4| is negligible relative to the OTHER coefficients' own scale, the
// polynomial is treated as whatever LOWER degree is actually present (recursing down through
// cubic/quadratic/linear/none) instead of dividing by ~0 -- this is what correctly handles the
// bN~=0 case above.
//
// Returns the number of real roots found (0-4, with multiplicity -- a genuine double root
// legitimately appears twice). Not sorted; roots[0..count) are written, the rest of roots[] is
// left untouched.
int solveQuarticReal(real_t a4, real_t a3, real_t a2, real_t a1, real_t a0, real_t roots[4]);

}  // namespace cardillo::solver
