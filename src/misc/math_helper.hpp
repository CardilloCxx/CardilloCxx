#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "types.hpp"

namespace cardillo {

    class MathHelper {
    public:
        static Quaternion4r alignQuaternionTo(const Quaternion4r& q_in, const Quaternion4r& q_ref) {
            Quaternion4r q = q_in;
            if (!q.coeffs().allFinite()) return q_ref;
            q.normalize();
            if (q_ref.dot(q) < (real_t)0) q.coeffs() = -q.coeffs();
            return q;
        }
    };

} // namespace cardillo
