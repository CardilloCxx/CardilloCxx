#pragma once

#include <cmath>
#include "misc/types.hpp"

namespace cardillo {
namespace collision {

// Build an orthonormal basis (t1, t2, n) from a unit-length normal n
// Uses a numerically stable, branchless variant based on Frisvad (2012).
// https://graphics.pixar.com/library/OrthonormalB/paper.pdf
inline void tangentFrameFromNormal(const Vector3r& n_in, Vector3r& t1, Vector3r& t2) {
    Vector3r n = n_in;
    real_t len = n.norm();
    if (len == (real_t)0) {
        // Degenerate input: default to world axes
        t1 = Vector3r::UnitX();
        t2 = Vector3r::UnitY();
        return;
    }
    n /= len;

    if (n.z() < (real_t)-0.9999999) {
        // Near south pole: choose fixed basis
        t1 = Vector3r(0, -1, 0);
        t2 = Vector3r(-1, 0, 0);
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
