#pragma once

#include <entt/entt.hpp>
#include <limits>
#include "../physics/api/physics_types.hpp"
#include "../physics/ecs_types.hpp"

namespace cardillo::RigidBody {

using RigidState = cardillo::physics::RigidState;

namespace detail {

inline RigidState readStateFromComponents(const entt::registry& reg, entt::entity e) {
    RigidState state = RigidState::inertial();
    if (!reg.valid(e)) return state;

    if (const auto* c = reg.try_get<C_Position3>(e)) state.position = c->value;
    if (const auto* c = reg.try_get<C_LinearVelocity3>(e)) state.linearVelocity = c->value;
    if (const auto* c = reg.try_get<C_AngularVelocity3>(e)) state.angularVelocity = c->value;
    if (const auto* c = reg.try_get<C_Orientation>(e)) state.orientation = c->value;
    state.rotation = state.orientation.toRotationMatrix();
    return state;
}

}  // namespace detail

inline RigidState getState(const entt::registry& reg, entt::entity e) {
    if (e == entt::null) return RigidState::inertial();
    if (const auto* cached = reg.try_get<RigidState>(e)) {
        return *cached;
    }
    return detail::readStateFromComponents(reg, e);
}

inline RigidState& updateState(entt::registry& reg, entt::entity e) {
    if (e == entt::null) {
        static RigidState inertial = RigidState::inertial();
        return inertial;
    }

    RigidState state = detail::readStateFromComponents(reg, e);
    return reg.emplace_or_replace<RigidState>(e, state);
}

}  // namespace cardillo::RigidBody