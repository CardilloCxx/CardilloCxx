#pragma once

#include <entt/entt.hpp>
#include "../physics/world.hpp"
#include "misc/types.hpp"

namespace cardillo {
namespace collision {

// Precomputed body transform data for fast world<->body conversions
struct BodyXform {
    Vector3r center{Vector3r::Zero()};   // world position of body origin
    Matrix33r R{Matrix33r::Identity()};  // rotation matrix: body -> world

    // Build from ECS components if available
    static BodyXform fromEcs(const entt::registry& reg, entt::entity e) {
        BodyXform X;
        if (reg.any_of<cardillo::C_Position3>(e)) X.center = reg.get<cardillo::C_Position3>(e).value;
        if (reg.any_of<cardillo::C_Orientation>(e)) X.R = reg.get<cardillo::C_Orientation>(e).value.toRotationMatrix();
        return X;
    }

    // World -> body
    Vector3r worldPointToBody(const Vector3r& pW) const { return R.transpose() * (pW - center); }
    Vector3r worldVecToBody(const Vector3r& vW) const { return R.transpose() * vW; }

    // Body -> world (not used here but provided for completeness)
    Vector3r bodyPointToWorld(const Vector3r& pB) const { return center + R * pB; }
    Vector3r bodyVecToWorld(const Vector3r& vB) const { return R * vB; }
};

}  // namespace collision
}  // namespace cardillo
