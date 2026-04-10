#include "physics_types.hpp"

#include "../world.hpp"

namespace cardillo::physics {

RigidState::RigidState(const Vector3r& p_local, const Vector3r& v_local, const Quaternion4r& q_local, const Vector3r& w_local, entt::entity refEntity, entt::registry& reg) {
    if (refEntity != entt::null && reg.all_of<cardillo::C_Position3, cardillo::C_Orientation, cardillo::C_LinearVelocity3, cardillo::C_AngularVelocity3>(refEntity)) {
        const auto& r_ORef = reg.get<cardillo::C_Position3>(refEntity).value;
        const auto& q_Ref = reg.get<cardillo::C_Orientation>(refEntity).value;
        const auto& v_Ref = reg.get<cardillo::C_LinearVelocity3>(refEntity).value;
        const auto& w_Ref = reg.get<cardillo::C_AngularVelocity3>(refEntity).value;

        const Matrix33r A_Ref = q_Ref.toRotationMatrix();
        position = r_ORef + A_Ref * p_local;
        orientation = q_Ref * q_local;

        const Vector3r r_rel_world = A_Ref * p_local;
        linearVelocity = v_Ref + w_Ref.cross(r_rel_world) + A_Ref * v_local;
        angularVelocity = w_Ref + A_Ref * w_local;
        return;
    }

    position = p_local;
    orientation = q_local;
    linearVelocity = v_local;
    angularVelocity = w_local;
}

}  // namespace cardillo::physics
