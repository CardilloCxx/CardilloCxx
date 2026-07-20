#pragma once

#include <entt/entt.hpp>
#include <string>
#include <utility>
#include <vector>
#include "../../misc/math_helper.hpp"
#include "beam_factory.hpp"
#include "rigid_body_factory.hpp"

#include "../world.hpp"

namespace cardillo {
namespace physics {

class BodyFactory {
   public:
    static entt::entity addRigidBody(World& sys, const physics::RigidShape& shape, const physics::RigidState& state, const physics::RigidProps& props) {
        return RigidBodyFactory::create(sys, shape, state, props);
    }

    static entt::entity addStaticBody(World& sys, const physics::RigidShape& shape, const physics::RigidState& state) {
        physics::RigidProps props(0);
        return RigidBodyFactory::create(sys, shape, state, props);
    }

    static std::pair<entt::entity, entt::entity> createBeam(World& sys, const misc::SplinePattern& spline, const physics::BeamCrossSection& section, const physics::BeamSpringParams& springs,
                                                            const physics::RigidState& stateDefaults, const physics::RigidProps& propsDefaults, size_t segments,
                                                            collision::CollisionCoal* collision_mgr = nullptr) {
        return BeamFactory::createBeam(sys, spline, section, springs, stateDefaults, propsDefaults, segments, collision_mgr);
    }

    static std::pair<entt::entity, entt::entity> createBeams(World& sys, const std::vector<const misc::SplinePattern*>& splines, const physics::BeamCrossSection& section,
                                                             const physics::BeamSpringParams& springs, const physics::RigidState& stateDefaults, const physics::RigidProps& propsDefaults,
                                                             size_t segments, collision::CollisionCoal* collision_mgr = nullptr) {
        return BeamFactory::createBeams(sys, splines, section, springs, stateDefaults, propsDefaults, segments, collision_mgr);
    }

    static entt::entity addPointMass(World& sys, real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius = (real_t)0.05);

    static std::vector<entt::entity> addSoftBody(World& sys, const std::string& objPath, real_t stiffness, real_t damping, const Vector3r& position = Vector3r::Zero(),
                                                 const Quaternion4r& orientation = Quaternion4r::Identity(), const Vector3r& linearVelocity = Vector3r::Zero(),
                                                 const Vector3r& angularVelocity = Vector3r::Zero(), real_t totalMass = (real_t)0.0, real_t nodeRadius = (real_t)0.02,
                                                 collision::CollisionCoal* collision_mgr = nullptr);
};

}  // namespace physics
}  // namespace cardillo
