#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include "../../misc/types.hpp"

namespace cardillo::solver::termination {

inline real_t velocityRms(const VectorXr& v, const VectorXr& v_prev, real_t tol_abs, real_t tol_rel) {
    if (v.size() == 0) return 0;
    const auto scale = v.array().abs().max(v_prev.array().abs()) * tol_rel + tol_abs;
    return 1.0 / std::sqrt((real_t)v.size()) * ((v - v_prev).array() / scale).matrix().norm();
}

// stepSize = per-row diagonal of R (not R^-1). Rows: frictionless normals first, then [N,T1,T2] triples.
// Returns NaN if anything is non-finite (caller maps to Diverged). Requires tol_abs > 0.
inline real_t naturalMapError(const VectorXr& lamOld, const VectorXr& lamNew, const VectorXr& stepSize,
                              const VectorXr& xi, const VectorXr& xiPrev, int frictionlessRows,
                              real_t tol_abs, real_t tol_rel) {
    real_t err = 0;
    bool finite = true;
    auto upd = [&](real_t r, real_t s) {
        const real_t e = std::abs(r) / (tol_abs + tol_rel * s);
        finite &= std::isfinite(e);
        err = std::max(err, e);
    };
    auto res = [&](int i) { return (lamOld[i] - lamNew[i]) / stepSize[i]; };
    auto vel = [&](int i) { return std::max(std::abs(xi[i]), std::abs(xiPrev[i])); };

    for (int i = 0; i < frictionlessRows; ++i) upd(res(i), vel(i));
    for (int i = frictionlessRows; i + 3 <= lamOld.size(); i += 3) {
        upd(res(i), vel(i));  // normal and tangent scaled separately
        upd(std::hypot(res(i + 1), res(i + 2)),
            std::max(std::hypot(xi[i + 1], xi[i + 2]), std::hypot(xiPrev[i + 1], xiPrev[i + 2])));
    }
    return finite ? err : std::numeric_limits<real_t>::quiet_NaN();
}

}  // namespace cardillo::solver::termination