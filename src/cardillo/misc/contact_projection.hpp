#pragma once

#include <cmath>
#include <limits>
#include "types.hpp"

namespace cardillo::misc {

// Contact laws in the negative variant of the projection equations (theory manual,
// eq:contact_law_prox / eq:ac_projection), shared by every projection-based solver
// (ProjectedJacobi, ProjectedGaussSeidel, Condensed incl. the local AC-Newton solve):
//
//   lambda_N = -proj_{R_0^-}(y_N),                 y_N = rho_N * xi_N - lambda_N
//   lambda_T = -proj_{D(mu * lambda_N)}(y_T),      y_T = rho_T * xi_T - lambda_T
//
// with xi the contact velocity, lambda the canonical contact percussion of `Contact::last_impulse`
// (lambda_N >= 0) and D(r) the disk of radius r. The new lambda_N defines
// the radius of the friction disk. Header-only and templated on the vector expression so that calls
// on fixed-size views inline completely -- no temporaries, no heap.

// proj_{R_0^-}(x)
inline real_t projNonPositive(real_t x) { return std::min(x, (real_t)0); }

// Frictionless contact: lambda_N = -proj_{R_0^-}(y_N).
inline real_t proxNormal(real_t yN) { return -projNonPositive(yN); }

// Frictional contact (3 rows): lambda = (-proj_{R_0^-}(y_N), -proj_{D(mu * lambda_N)}(y_T)).
// `y` may alias `lambda` (it is copied first), so `proxCoulomb(step - lam, mu, lam)` is fine.
template <typename In, typename Out>
inline void proxCoulomb(const Eigen::MatrixBase<In>& y_in, real_t mu, Out&& lambda) {
    const Vector3r y = y_in;
    lambda[0] = proxNormal(y[0]);

    // -proj_D(y_T): y_T scaled onto the disk of radius mu * lambda_N, then negated.
    const real_t yT_norm_sq = y[1] * y[1] + y[2] * y[2];
    const real_t s = (yT_norm_sq > (real_t)0) ? std::min((real_t)1, (mu * lambda[0]) / std::sqrt(yT_norm_sq)) : (real_t)1;
    lambda[1] = -s * y[1];
    lambda[2] = -s * y[2];
}

// Stopping test of the projection solvers. The increment of one projection step above, divided by
// the step rho of the contact,
//
//   r = (lambda - lambda^+) / rho,      lambda^+ the result of proxNormal / proxCoulomb,
//
// is the natural map of the contact law at the input (lambda, xi) of the step: a velocity that is
// zero iff (lambda, xi) satisfies the contact law (both tangential rows need the same rho), for a
// sticking contact the velocity error itself. naturalMapError2() returns e^2 of
//
//   e = max(|r_N| / (tolAbs + tolRel |xi_N|), ||r_T|| / (tolAbs + tolRel ||xi_T||)),
//
// e <= 1: within tolerance. Normal and tangential parts are scaled separately, so that a normal
// error is not hidden behind a large slip. tolAbs > 0. A NaN error is returned as infinity, which
// survives the maximum over the contacts. The frictional variant takes the rows (n, t1, t2) as
// scalars: building them as Vector3r expressions in the projection loop costs three times as much.
inline real_t naturalMapError2(real_t r, real_t xi, real_t tolAbs, real_t tolRel) {
    const real_t d = tolAbs + tolRel * std::abs(xi);
    const real_t e2 = r * r / (d * d);
    return std::isnan(e2) ? std::numeric_limits<real_t>::infinity() : e2;
}
inline real_t naturalMapError2(real_t rN, real_t rT1, real_t rT2, real_t xiN, real_t xiT1, real_t xiT2, real_t tolAbs, real_t tolRel) {
    const real_t dN = tolAbs + tolRel * std::abs(xiN);
    const real_t dT = tolAbs + tolRel * std::sqrt(xiT1 * xiT1 + xiT2 * xiT2);
    const real_t eN2 = rN * rN / (dN * dN);
    const real_t eT2 = (rT1 * rT1 + rT2 * rT2) / (dT * dT);
    return std::isnan(eN2 + eT2) ? std::numeric_limits<real_t>::infinity() : std::max(eN2, eT2);
}

// Euclidean projection onto the Coulomb cone and the De Saxce correction, for one frictional
// contact in the row order (n, t1, t2), used by the natural map of the final check
// (solver/natural_map.hpp):
//
//   K_mu  = { x : ||x_T|| <= mu x_N }          (percussions)
//   K_mu* = { v : mu ||v_T|| <= v_N }          (corrected relative velocities)
//   Gamma(xi) = (mu ||xi_T||, 0, 0)            (De Saxce: pi in K_mu perp xi + Gamma(xi) in K_mu*)
//
// The polar cone of K_mu is -K_mu*, so x = proj_K(x) + proj_{-K*}(x) with orthogonal parts
// (Moreau decomposition). Both functions are positively homogeneous of degree one, hence invariant
// under one positive scaling factor per contact.

// proj_{K_mu}(x). The polar test comes first: at mu = 0 it catches every point with x_N <= 0.
inline Vector3r projectCoulombCone(const Vector3r& x, real_t mu) {
    const real_t tNorm = std::sqrt(x[1] * x[1] + x[2] * x[2]);
    if (mu * tNorm <= -x[0]) return Vector3r::Zero();
    if (tNorm <= mu * x[0]) return x;
    // Boundary point (r, mu r e_T) closest to x, r = (x_N + mu ||x_T||) / (1 + mu^2); tNorm > 0 here.
    const real_t r = (x[0] + mu * tNorm) / ((real_t)1 + mu * mu);
    const real_t scaleT = mu * r / tNorm;
    return Vector3r(r, scaleT * x[1], scaleT * x[2]);
}

// Normal component of the De Saxce correction Gamma(xi) = (mu ||xi_T||, 0, 0).
inline real_t deSaxce(const Vector3r& xi, real_t mu) {
    return mu * std::sqrt(xi[1] * xi[1] + xi[2] * xi[2]);
}

}  // namespace cardillo::misc
