#pragma once

#include <entt/entt.hpp>
#include <limits>
#include "../physics/api/physics_types.hpp"
#include "../physics/ecs_types.hpp"

/**
 * @brief Read-only accessors for a rigid body's kinematic state (see cardillo::physics::RigidState
 * for the frame convention: world-frame position/orientation/linear velocity, body-frame angular
 * velocity).
 */
namespace cardillo::RigidBody {

using RigidState = cardillo::physics::RigidState;

namespace detail {

/// Reconstructs a RigidState from its individual ECS components; used when no cached RigidState exists.
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

/**
 * @brief Returns the current kinematic state of entity @p e, preferring the cached RigidState
 * component if present (see updateState()), otherwise reconstructing it from raw components.
 * @return `RigidState::inertial()` (identity pose, zero velocity) for entt::null or an entity with
 * no position/orientation components.
 */
inline RigidState getState(const entt::registry& reg, entt::entity e) {
    if (e == entt::null) return RigidState::inertial();
    if (const auto* cached = reg.try_get<RigidState>(e)) {
        return *cached;
    }
    return detail::readStateFromComponents(reg, e);
}

/// An entity is static if it lacks a dynamics index/physics-object tag (fixed bodies, or invalid entities).
inline bool isStatic(const entt::registry& reg, entt::entity e) {
    if (!reg.valid(e)) return true;  // Non-existent entities are considered static
    return !(reg.any_of<C_BodyIndex, C_PhysicsObject>(e));
}

/**
 * @brief Recomputes and caches entity @p e's RigidState component from its raw position/velocity
 * components, replacing any previously cached value.
 * @warning For `e == entt::null`, this returns a reference to a shared, mutable, function-local
 * static fallback object rather than a per-call value; callers must not write through it.
 */
inline RigidState& updateState(entt::registry& reg, entt::entity e) {
    if (e == entt::null) {
        static RigidState inertial = RigidState::inertial();
        return inertial;
    }

    RigidState state = detail::readStateFromComponents(reg, e);
    return reg.emplace_or_replace<RigidState>(e, state);
}

}  // namespace cardillo::RigidBody