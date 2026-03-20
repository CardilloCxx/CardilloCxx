#pragma once

#include <string>
#include <vector>
#include <utility>
#include <entt/entt.hpp>

#include "physics_system.hpp"

namespace cardillo {
namespace physics {

class BodyFactory {
public:
    static entt::entity addRigidBody(PhysicsSystem& sys,
                                     const PhysicsSystem::RigidShape& shape,
                                     const PhysicsSystem::RigidState& state,
                                     const PhysicsSystem::RigidProps& props);

    static entt::entity addStaticBody(PhysicsSystem& sys,
                                      const PhysicsSystem::RigidShape& shape,
                                      const PhysicsSystem::RigidState& state);

    static index_t addPointMass(PhysicsSystem& sys,
                                real_t mass,
                                const Vector3r& x0,
                                const Vector3r& v0,
                                real_t radius = (real_t)0.05);

    static index_t addObstacleHeightField(PhysicsSystem& sys,
                                          const Vector3r& position,
                                          const Quaternion4r& orientation,
                                          const std::string& exrPath,
                                          real_t x_dim,
                                          real_t y_dim,
                                          real_t z_scale = (real_t)1.0,
                                          real_t min_height = (real_t)0.0);

    static std::vector<entt::entity> addSoftBody(PhysicsSystem& sys,
                                                  const std::string& objPath,
                                                  real_t stiffness,
                                                  real_t damping,
                                                  const Vector3r& position = Vector3r::Zero(),
                                                  const Quaternion4r& orientation = Quaternion4r::Identity(),
                                                  const Vector3r& linearVelocity = Vector3r::Zero(),
                                                  const Vector3r& angularVelocity = Vector3r::Zero(),
                                                  real_t totalMass = (real_t)0.0);

    static std::pair<entt::entity, entt::entity> createBeam(PhysicsSystem& sys,
                                                             const misc::SplinePattern& spline,
                                                             const PhysicsSystem::BeamCrossSection& section,
                                                             const PhysicsSystem::BeamSpringParams& springs,
                                                             const PhysicsSystem::RigidState& stateDefaults,
                                                             const PhysicsSystem::RigidProps& propsDefaults,
                                                             size_t segments);

    static std::pair<entt::entity, entt::entity> createBeams(PhysicsSystem& sys,
                                                              const std::vector<const misc::SplinePattern*>& splines,
                                                              const PhysicsSystem::BeamCrossSection& section,
                                                              const PhysicsSystem::BeamSpringParams& springs,
                                                              const PhysicsSystem::RigidState& stateDefaults,
                                                              const PhysicsSystem::RigidProps& propsDefaults,
                                                              size_t segmentsPerSpline);
};

} // namespace physics
} // namespace cardillo
