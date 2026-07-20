#include "physics_types.hpp"

#include "../../rigid_body/transformations.hpp"
#include "../world.hpp"

namespace cardillo::physics {

RigidState::RigidState(const Vector3r& p_local, const Vector3r& v_local, const Quaternion4r& q_local, const Vector3r& w_local, entt::entity refEntity, entt::registry& reg) {
    const RigidBody::RigidState refState = RigidBody::getState(reg, refEntity);
    const RigidBody::RigidState inertial = RigidBody::RigidState::inertial();

    RigidState localState;
    localState.position = p_local;
    localState.orientation = q_local;
    localState.rotation = q_local.toRotationMatrix();
    localState.linearVelocity = v_local;
    localState.angularVelocity = w_local;

    *this = transform::rigidState(localState, refState, inertial);
}

}  // namespace cardillo::physics
