#pragma once

#include <cmath>
#include <limits>
#include "misc/types.hpp"

namespace cardillo {
namespace collision {

/**
 * @brief Builds a right-handed orthonormal basis (t1, t2, n) from a contact normal.
 *
 * Uses a numerically stable variant of Frisvad's construction
 * (https://graphics.pixar.com/library/OrthonormalB/paper.pdf), 
 * (https://backend.orbit.dtu.dk/ws/files/126824972/onb_frisvad_jgt2012_v2.pdf), with an additional
 * re-orthogonalization step near the south pole (n.z() -> -1) to keep (t1, t2) exactly
 * orthogonal to `n_in` even when it is only approximately antiparallel to +z.
 *
 * @param n_in Contact normal; need not be pre-normalized (any nonzero length is accepted).
 * @param[out] t1 First tangent direction, orthogonal to the normalized normal.
 * @param[out] t2 Second tangent direction, such that (t1, t2, n) is right-handed.
 */
inline void tangentFrameFromNormal(const Vector3r& n_in, Vector3r& t1, Vector3r& t2) {
    Vector3r n = n_in;
    real_t len = n.norm();
    if (len <= std::numeric_limits<real_t>::epsilon()) {
        // Degenerate input: default to world axes
        t1 = Vector3r::UnitX();
        t2 = Vector3r::UnitY();
        return;
    }
    n /= len;

    if (n.z() < (real_t)-0.9999999) {
        // Near south pole: start from a fixed basis, then re-orthogonalize against the actual
        // (near-pole but not exact) normal so (t1, t2, n) stays right-handed and orthonormal.
        t1 = Vector3r(0, -1, 0);
        t1 -= n * n.dot(t1);
        t1.normalize();
        t2 = n.cross(t1);
        return;
    }
    real_t a = (real_t)1 / ((real_t)1 + n.z());
    real_t b = -n.x() * n.y() * a;
    t1 = Vector3r((real_t)1 - n.x() * n.x() * a, b, -n.x());
    t2 = Vector3r(b, (real_t)1 - n.y() * n.y() * a, -n.y());
    // Optional: normalize to guard against accumulated error
    t1.normalize();
    t2 = n.cross(t1);  // ensure right-handed ONB
}

}  // namespace collision
}  // namespace cardillo
