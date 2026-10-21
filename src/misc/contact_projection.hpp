#pragma once

#include <algorithm>
#include <cmath>
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
// (see solver/warmstart.hpp; lambda_N >= 0) and D(r) the disk of radius r. The new lambda_N defines
// the radius of the friction disk. Header-only and templated on the vector expression so that calls
// on fixed-size views inline completely -- no temporaries, no heap.

// proj_{R_0^-}(x)
inline real_t projNonPositive(real_t x) { return std::min((real_t)0, x); }

// Frictionless contact: lambda_N = -proj_{R_0^-}(y_N).
inline real_t proxNormal(real_t yN) { return -projNonPositive(yN); }

// Frictional contact (3 rows): lambda = (-proj_{R_0^-}(y_N), -proj_{D(mu * lambda_N)}(y_T)).
// `y` may alias `lambda` (it is copied first), so `proxCoulomb(step - lam, mu, lam)` is fine.
template <typename In, typename Out>
inline void proxCoulomb(const Eigen::MatrixBase<In>& y_in, real_t mu, Out&& lambda) {
    const Vector3r y = y_in;
    lambda[0] = -projNonPositive(y[0]);

    // -proj_D(y_T): y_T scaled onto the disk of radius mu * lambda_N, then negated.
    const real_t yT_norm_sq = y[1] * y[1] + y[2] * y[2];
    const real_t s = (yT_norm_sq > (real_t)0) ? std::min((real_t)1, (mu * lambda[0]) / std::sqrt(yT_norm_sq)) : (real_t)1;
    lambda[1] = -s * y[1];
    lambda[2] = -s * y[2];
}

}  // namespace cardillo::misc
