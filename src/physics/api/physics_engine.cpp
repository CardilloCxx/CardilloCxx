#include "physics_engine.hpp"

#include "../world.hpp"
#include "../assets/body_factory.hpp"
#include "../constraints/constraint_factory.hpp"

namespace cardillo {
namespace physics {

// World now aliases public physics types; no conversion helpers required.

// Trivial constructor and accessors are inlined in the header.

entt::entity PhysicsEngine::addRigidBody(const RigidShape& shape,
                                         const RigidState& state,
                                         const RigidProps& props) {
    return BodyFactory::addRigidBody(m_world, shape, state, props);
}

entt::entity PhysicsEngine::addStaticBody(const RigidShape& shape,
                                          const RigidState& state) {
    return BodyFactory::addStaticBody(m_world, shape, state);
}

index_t PhysicsEngine::addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius) {
    return BodyFactory::addPointMass(m_world, mass, x0, v0, radius);
}

index_t PhysicsEngine::addObstacleHeightField(const Vector3r& position,
                                             const Quaternion4r& orientation,
                                             const std::string& exrPath,
                                             real_t x_dim,
                                             real_t y_dim,
                                             real_t z_scale,
                                             real_t min_height) {
    return BodyFactory::addObstacleHeightField(m_world, position, orientation, exrPath, x_dim, y_dim, z_scale, min_height);
}

std::vector<entt::entity> PhysicsEngine::addSoftBody(const std::string& objPath,
                                                    real_t stiffness,
                                                    real_t damping,
                                                    const Vector3r& position,
                                                    const Quaternion4r& orientation,
                                                    const Vector3r& linearVelocity,
                                                    const Vector3r& angularVelocity,
                                                    real_t totalMass) {
    return BodyFactory::addSoftBody(m_world, objPath, stiffness, damping, position, orientation, linearVelocity, angularVelocity, totalMass);
}

std::pair<entt::entity, entt::entity> PhysicsEngine::createBeam(const cardillo::misc::SplinePattern& spline,
                                                                const BeamCrossSection& section,
                                                                const BeamSpringParams& springs,
                                                                const RigidState& stateDefaults,
                                                                const RigidProps& propsDefaults,
                                                                size_t segments) {
    return BodyFactory::createBeam(m_world, spline, section, springs, stateDefaults, propsDefaults, segments);
}

std::pair<entt::entity, entt::entity> PhysicsEngine::createBeams(const std::vector<const cardillo::misc::SplinePattern*>& splines,
                                                                 const BeamCrossSection& section,
                                                                 const BeamSpringParams& springs,
                                                                 const RigidState& stateDefaults,
                                                                 const RigidProps& propsDefaults,
                                                                 size_t segmentsPerSpline) {
    return BodyFactory::createBeams(m_world, splines, section, springs, stateDefaults, propsDefaults, segmentsPerSpline);
}

size_t PhysicsEngine::addLinearDistanceConstraint(entt::entity a, entt::entity b, const Vector3r& rA_local, const Vector3r& rB_local, real_t stiffness, real_t damping) {
    return ConstraintFactory::addLinearDistanceConstraint(m_world, a, b, rA_local, rB_local, stiffness, damping);
}

size_t PhysicsEngine::addRigidConstraint(entt::entity a, entt::entity b) {
    return ConstraintFactory::addRigidConstraint(m_world, a, b);
}

size_t PhysicsEngine::addTranslationRotationConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans, const Vector3r& D_trans, const Vector3r& K_rot, const Vector3r& D_rot) {
    return ConstraintFactory::addTranslationRotationConstraint(m_world, a, b, frame, K_trans, D_trans, K_rot, D_rot);
}

size_t PhysicsEngine::addTranslationalConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans, const Vector3r& D_trans) {
    return ConstraintFactory::addTranslationalConstraint(m_world, a, b, frame, K_trans, D_trans);
}

size_t PhysicsEngine::addRotationConstraint(entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_rot, const Vector3r& D_rot) {
    return ConstraintFactory::addRotationConstraint(m_world, a, b, frame, K_rot, D_rot);
}

size_t PhysicsEngine::addHingeConstraint(entt::entity a, entt::entity b, const JointFrame& frame, real_t K_axis, real_t D_axis, const Vector3r& K_trans, const Vector3r& D_trans) {
    return ConstraintFactory::addHingeConstraint(m_world, a, b, frame, K_axis, D_axis, K_trans, D_trans);
}

size_t PhysicsEngine::addBeamConstraint(entt::entity a, entt::entity b, const BeamSpringParams& springs, const BeamCrossSection& section) {
    return ConstraintFactory::addBeamConstraint(m_world, a, b, springs, section);
}

} // namespace physics
} // namespace cardillo
