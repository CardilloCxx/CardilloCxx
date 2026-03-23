#include "physics_engine.hpp"

#include "../world.hpp"
#include "../assets/body_factory.hpp"
#include "../constraints/constraint_factory.hpp"

namespace cardillo {
namespace physics {

// Helper converters between api types and World types
static cardillo::World::RigidState toWorld(const RigidState& s) {
    cardillo::World::RigidState w;
    w.position = s.position;
    w.orientation = s.orientation;
    w.linearVelocity = s.linearVelocity;
    w.angularVelocity = s.angularVelocity;
    return w;
}

static cardillo::World::RigidProps toWorld(const RigidProps& p) {
    cardillo::World::RigidProps w;
    w.mass = p.mass;
    w.density = p.density;
    w.friction = p.friction;
    w.collidable = p.collidable;
    w.visual = p.visual;
    return w;
}

static cardillo::World::RigidShape toWorld(const RigidShape& s) {
    return std::visit([](auto&& arg) -> cardillo::World::RigidShape {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, CubeShape>) return cardillo::World::RigidShape{cardillo::World::CubeShape{arg.halfExtents}};
        else if constexpr (std::is_same_v<T, PlaneShape>) return cardillo::World::RigidShape{cardillo::World::PlaneShape{arg.normal, arg.up, arg.sizeX, arg.sizeY}};
        else if constexpr (std::is_same_v<T, CapsuleShape>) return cardillo::World::RigidShape{cardillo::World::CapsuleShape{arg.radius, arg.halfLength}};
        else if constexpr (std::is_same_v<T, CylinderShape>) return cardillo::World::RigidShape{cardillo::World::CylinderShape{arg.radius, arg.halfLength}};
        else if constexpr (std::is_same_v<T, ConeShape>) return cardillo::World::RigidShape{cardillo::World::ConeShape{arg.radius, arg.height}};
        else if constexpr (std::is_same_v<T, SphereShape>) return cardillo::World::RigidShape{cardillo::World::SphereShape{arg.radius}};
        else if constexpr (std::is_same_v<T, MeshShape>) return cardillo::World::RigidShape{cardillo::World::MeshShape{arg.path, arg.scale, arg.use_bbox_collider, arg.show_collider}};
        else return cardillo::World::RigidShape{};
    }, s);
}

static cardillo::World::BeamCrossSection toWorld(const BeamCrossSection& s) {
    return cardillo::World::BeamCrossSection{s.width, s.height, static_cast<cardillo::World::BeamBodyType>(s.type)};
}

static cardillo::World::BeamSpringParams toWorld(const BeamSpringParams& p) {
    cardillo::World::BeamSpringParams w;
    w.E = p.E; w.nu = p.nu; w.scaleKe = p.scaleKe; w.scaleKf = p.scaleKf; w.Ke_direct = p.Ke_direct; w.Kf_direct = p.Kf_direct; w.dampingFactor = p.dampingFactor;
    return w;
}

// Trivial constructor and accessors are inlined in the header.

entt::entity PhysicsEngine::addRigidBody(const RigidShape& shape,
                                         const RigidState& state,
                                         const RigidProps& props) {
    return BodyFactory::addRigidBody(m_world, toWorld(shape), toWorld(state), toWorld(props));
}

entt::entity PhysicsEngine::addStaticBody(const RigidShape& shape,
                                          const RigidState& state) {
    return BodyFactory::addStaticBody(m_world, toWorld(shape), toWorld(state));
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
    return BodyFactory::createBeam(m_world, spline, toWorld(section), toWorld(springs), toWorld(stateDefaults), toWorld(propsDefaults), segments);
}

std::pair<entt::entity, entt::entity> PhysicsEngine::createBeams(const std::vector<const cardillo::misc::SplinePattern*>& splines,
                                                                 const BeamCrossSection& section,
                                                                 const BeamSpringParams& springs,
                                                                 const RigidState& stateDefaults,
                                                                 const RigidProps& propsDefaults,
                                                                 size_t segmentsPerSpline) {
    return BodyFactory::createBeams(m_world, splines, toWorld(section), toWorld(springs), toWorld(stateDefaults), toWorld(propsDefaults), segmentsPerSpline);
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
    return ConstraintFactory::addBeamConstraint(m_world, a, b, toWorld(springs), toWorld(section));
}

entt::registry& PhysicsEngine::ecs() { return m_world.ecs(); }
const entt::registry& PhysicsEngine::ecs() const { return m_world.ecs(); }

void PhysicsEngine::disableCollisionBetween(entt::entity a, entt::entity b) { m_world.disableCollisionBetween(a,b); }
void PhysicsEngine::makeStatic(entt::entity e) { m_world.makeStatic(e); }

MatrixXXr PhysicsEngine::getMass(entt::entity e) const { return m_world.getMass(e); }
Vector3r PhysicsEngine::getInertiaDiag(entt::entity e) const { return m_world.getInertiaDiag(e); }
VectorXr PhysicsEngine::getPosition(entt::entity e) const { return m_world.getPosition(e); }
real_t PhysicsEngine::getKineticEnergy(entt::entity e) const { return m_world.getKineticEnergy(e); }

void PhysicsEngine::applyForce(entt::entity e, const Vector3r& f, const Vector3r& tau) { m_world.applyForce(e, f, tau); }
void PhysicsEngine::applyInertialTorque(entt::entity e, const Vector3r& tau) { m_world.applyInertialTorque(e, tau); }

void PhysicsEngine::setPosition(entt::entity e, const Vector3r& p) { m_world.setPosition(e, p); }
void PhysicsEngine::setOrientation(entt::entity e, const Quaternion4r& q) { m_world.setOrientation(e, q); }
void PhysicsEngine::setLinearVelocity(entt::entity e, const Vector3r& v) { m_world.setLinearVelocity(e, v); }
void PhysicsEngine::setAngularVelocity(entt::entity e, const Vector3r& w) { m_world.setAngularVelocity(e, w); }
void PhysicsEngine::setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w) { m_world.setVelocityByForce(e, v, w); }

void PhysicsEngine::setGravity(const Vector3r& g) { m_world.setGravity(g); }
const Vector3r& PhysicsEngine::gravity() const { return m_world.gravity(); }

void PhysicsEngine::markStructureDirty() { m_world.markStructureDirty(); }
void PhysicsEngine::track(entt::entity e, const std::string& name) { m_world.track(e, name); }
void PhysicsEngine::writeTrackedStateToCsv(real_t t) { m_world.writeTrackedStateToCsv(t); }

} // namespace physics
} // namespace cardillo
