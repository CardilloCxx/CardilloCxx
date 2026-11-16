#pragma once

#include <cmath>
#include "../../misc/types.hpp"

namespace cardillo::solver::termination {

// RMS of the generalized-velocity increment, each component scaled by tol_abs + tol_rel * max(|v|, |v_prev|)
// (pj.criterion = velocity_rms); <= 1: converged. The alternative pj.criterion = natural_map is a
// residual of the contact law (misc::naturalMapError2).
inline real_t velocityRms(const VectorXr& v, const VectorXr& v_prev, real_t tol_abs, real_t tol_rel) {
    if (v.size() == 0) return 0;
    const auto scale = v.array().abs().max(v_prev.array().abs()) * tol_rel + tol_abs;
    return 1.0 / std::sqrt((real_t)v.size()) * ((v - v_prev).array() / scale).matrix().norm();
}

}  // namespace cardillo::solver::termination
