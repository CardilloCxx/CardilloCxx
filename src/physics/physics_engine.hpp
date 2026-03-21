#pragma once

#include <type_traits>

#include "physics_types.hpp"
#include "body_factory.hpp"
#include "constraint_factory.hpp"

namespace cardillo {
namespace physics {

class PhysicsEngine {
public:
    explicit PhysicsEngine(World& world) : m_world(world) {}

    World& world() { return m_world; }
    const World& world() const { return m_world; }

    World& system() { return m_world; }
    const World& system() const { return m_world; }

    entt::entity addRigidBody(const RigidShape& shape,
                              const RigidState& state,
                              const RigidProps& props) {
        return BodyFactory::addRigidBody(m_world, toWorldShape(shape), toWorldState(state), toWorldProps(props));
    }

    entt::entity addRigidBody(const World::RigidShape& shape,
                              const World::RigidState& state,
                              const World::RigidProps& props) {
        return BodyFactory::addRigidBody(m_world, shape, state, props);
    }

    entt::entity addStaticBody(const RigidShape& shape,
                               const RigidState& state) {
        return BodyFactory::addStaticBody(m_world, toWorldShape(shape), toWorldState(state));
    }

    entt::entity addStaticBody(const World::RigidShape& shape,
                               const World::RigidState& state) {
        return BodyFactory::addStaticBody(m_world, shape, state);
    }

    index_t addPointMass(real_t mass,
                         const Vector3r& x0,
                         const Vector3r& v0,
                         real_t radius = (real_t)0.05) {
        return BodyFactory::addPointMass(m_world, mass, x0, v0, radius);
    }

    index_t addObstacleHeightField(const Vector3r& position,
                                   const Quaternion4r& orientation,
                                   const std::string& exrPath,
                                   real_t x_dim,
                                   real_t y_dim,
                                   real_t z_scale = (real_t)1.0,
                                   real_t min_height = (real_t)0.0) {
        return BodyFactory::addObstacleHeightField(m_world, position, orientation, exrPath, x_dim, y_dim, z_scale, min_height);
    }

    std::vector<entt::entity> addSoftBody(const std::string& objPath,
                                          real_t stiffness,
                                          real_t damping,
                                          const Vector3r& position = Vector3r::Zero(),
                                          const Quaternion4r& orientation = Quaternion4r::Identity(),
                                          const Vector3r& linearVelocity = Vector3r::Zero(),
                                          const Vector3r& angularVelocity = Vector3r::Zero(),
                                          real_t totalMass = (real_t)0.0) {
        return BodyFactory::addSoftBody(m_world, objPath, stiffness, damping, position, orientation, linearVelocity, angularVelocity, totalMass);
    }

    std::pair<entt::entity, entt::entity> createBeam(const misc::SplinePattern& spline,
                                                      const BeamCrossSection& section,
                                                      const BeamSpringParams& springs,
                                                      const RigidState& stateDefaults,
                                                      const RigidProps& propsDefaults,
                                                      size_t segments) {
        return BodyFactory::createBeam(
            m_world,
            spline,
            toWorldSection(section),
            toWorldSprings(springs),
            toWorldState(stateDefaults),
            toWorldProps(propsDefaults),
            segments);
    }

    std::pair<entt::entity, entt::entity> createBeam(const misc::SplinePattern& spline,
                                                      const World::BeamCrossSection& section,
                                                      const World::BeamSpringParams& springs,
                                                      const World::RigidState& stateDefaults,
                                                      const World::RigidProps& propsDefaults,
                                                      size_t segments) {
        return BodyFactory::createBeam(m_world, spline, section, springs, stateDefaults, propsDefaults, segments);
    }

    std::pair<entt::entity, entt::entity> createBeams(const std::vector<const misc::SplinePattern*>& splines,
                                                       const BeamCrossSection& section,
                                                       const BeamSpringParams& springs,
                                                       const RigidState& stateDefaults,
                                                       const RigidProps& propsDefaults,
                                                       size_t segmentsPerSpline) {
        return BodyFactory::createBeams(
            m_world,
            splines,
            toWorldSection(section),
            toWorldSprings(springs),
            toWorldState(stateDefaults),
            toWorldProps(propsDefaults),
            segmentsPerSpline);
    }

    std::pair<entt::entity, entt::entity> createBeams(const std::vector<const misc::SplinePattern*>& splines,
                                                       const World::BeamCrossSection& section,
                                                       const World::BeamSpringParams& springs,
                                                       const World::RigidState& stateDefaults,
                                                       const World::RigidProps& propsDefaults,
                                                       size_t segmentsPerSpline) {
        return BodyFactory::createBeams(m_world, splines, section, springs, stateDefaults, propsDefaults, segmentsPerSpline);
    }

    size_t addLinearDistanceConstraint(entt::entity a,
                                       entt::entity b,
                                       const Vector3r& rA_local = Vector3r::Zero(),
                                       const Vector3r& rB_local = Vector3r::Zero(),
                                       real_t stiffness = std::numeric_limits<real_t>::infinity(),
                                       real_t damping = (real_t)0) {
        return ConstraintFactory::addLinearDistanceConstraint(m_world, a, b, rA_local, rB_local, stiffness, damping);
    }

    size_t addRigidConstraint(entt::entity a, entt::entity b) {
        return ConstraintFactory::addRigidConstraint(m_world, a, b);
    }

    size_t addTranslationRotationConstraint(entt::entity a,
                                            entt::entity b,
                                            const JointFrame& frame,
                                            const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                            const Vector3r& D_trans = Vector3r::Zero(),
                                            const Vector3r& K_rot = Vector3r::Zero(),
                                            const Vector3r& D_rot = Vector3r::Zero()) {
        return ConstraintFactory::addTranslationRotationConstraint(m_world, a, b, frame, K_trans, D_trans, K_rot, D_rot);
    }

    size_t addTranslationalConstraint(entt::entity a,
                                      entt::entity b,
                                      const JointFrame& frame,
                                      const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                      const Vector3r& D_trans = Vector3r::Zero()) {
        return ConstraintFactory::addTranslationalConstraint(m_world, a, b, frame, K_trans, D_trans);
    }

    size_t addRotationConstraint(entt::entity a,
                                 entt::entity b,
                                 const JointFrame& frame,
                                 const Vector3r& K_rot = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                 const Vector3r& D_rot = Vector3r::Zero()) {
        return ConstraintFactory::addRotationConstraint(m_world, a, b, frame, K_rot, D_rot);
    }

    size_t addHingeConstraint(entt::entity a,
                              entt::entity b,
                              const JointFrame& frame,
                              real_t K_axis = (real_t)0,
                              real_t D_axis = (real_t)0,
                              const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                              const Vector3r& D_trans = Vector3r::Zero()) {
        return ConstraintFactory::addHingeConstraint(m_world, a, b, frame, K_axis, D_axis, K_trans, D_trans);
    }

    size_t addBeamConstraint(entt::entity a,
                             entt::entity b,
                             const BeamSpringParams& springs,
                             const BeamCrossSection& section) {
        return ConstraintFactory::addBeamConstraint(m_world, a, b, toWorldSprings(springs), toWorldSection(section));
    }

    size_t addBeamConstraint(entt::entity a,
                             entt::entity b,
                             const World::BeamSpringParams& springs,
                             const World::BeamCrossSection& section) {
        return ConstraintFactory::addBeamConstraint(m_world, a, b, springs, section);
    }

private:
    static World::RigidState toWorldState(const RigidState& s) {
        World::RigidState out;
        out.position = s.position;
        out.orientation = s.orientation;
        out.linearVelocity = s.linearVelocity;
        out.angularVelocity = s.angularVelocity;
        return out;
    }

    static World::RigidProps toWorldProps(const RigidProps& p) {
        World::RigidProps out;
        out.mass = p.mass;
        out.density = p.density;
        out.friction = p.friction;
        out.collidable = p.collidable;
        out.visual = p.visual;
        return out;
    }

    static World::BeamCrossSection toWorldSection(const BeamCrossSection& s) {
        World::BeamBodyType type = World::BeamBodyType::Cube;
        switch (s.type) {
            case BeamBodyType::Cube: type = World::BeamBodyType::Cube; break;
            case BeamBodyType::Capsule: type = World::BeamBodyType::Capsule; break;
            case BeamBodyType::Cylinder: type = World::BeamBodyType::Cylinder; break;
        }
        return World::BeamCrossSection(s.width, s.height, type);
    }

    static World::BeamSpringParams toWorldSprings(const BeamSpringParams& s) {
        World::BeamSpringParams out;
        out.E = s.E;
        out.nu = s.nu;
        out.scaleKe = s.scaleKe;
        out.scaleKf = s.scaleKf;
        out.Ke_direct = s.Ke_direct;
        out.Kf_direct = s.Kf_direct;
        out.dampingFactor = s.dampingFactor;
        return out;
    }

    static World::RigidShape toWorldShape(const RigidShape& shape) {
        return std::visit(
            [](const auto& s) -> World::RigidShape {
                using T = std::decay_t<decltype(s)>;
                if constexpr (std::is_same_v<T, CubeShape>) {
                    return World::CubeShape{s.halfExtents};
                } else if constexpr (std::is_same_v<T, PlaneShape>) {
                    return World::PlaneShape{s.normal, s.up, s.sizeX, s.sizeY};
                } else if constexpr (std::is_same_v<T, CapsuleShape>) {
                    return World::CapsuleShape{s.radius, s.halfLength};
                } else if constexpr (std::is_same_v<T, CylinderShape>) {
                    return World::CylinderShape{s.radius, s.halfLength};
                } else if constexpr (std::is_same_v<T, ConeShape>) {
                    return World::ConeShape{s.radius, s.height};
                } else if constexpr (std::is_same_v<T, SphereShape>) {
                    return World::SphereShape{s.radius};
                } else {
                    return World::MeshShape{s.path, s.scale, s.use_bbox_collider, s.show_collider};
                }
            },
            shape);
    }

    World& m_world;
};

} // namespace physics
} // namespace cardillo
