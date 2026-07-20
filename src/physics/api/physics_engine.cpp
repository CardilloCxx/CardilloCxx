#include "physics_engine.hpp"

#include "../../collision/collision_coal.hpp"
#include "../../misc/timings/TimingManager.hpp"
#include "../pipeline/physics_pipeline.hpp"

namespace cardillo {
namespace physics {

PhysicsEngine::~PhysicsEngine() = default;

PhysicsEngine::PhysicsEngine() = default;

PhysicsEngine::PhysicsEngine(const config::Config& cfg) : PhysicsEngine() {
    initFromConfig(cfg);
}

void PhysicsEngine::initFromConfig(const config::Config& cfg) {
    m_cfg = cfg;
    if (!m_world) m_world = std::make_unique<World>(m_cfg);
    m_timings = std::make_unique<misc::TimingManager>();
    m_collision_mgr = std::make_unique<collision::CollisionCoal>(*m_world, m_timings.get(), m_cfg);
    m_pipeline = std::make_unique<physics::pipeline::PhysicsPipeline>(*m_world, m_cfg, m_collision_mgr.get(), m_timings.get());
}

void PhysicsEngine::step(real_t dt) {
    if (m_pipeline) m_pipeline->step(dt);
}

void PhysicsEngine::step() {
    if (m_pipeline) m_pipeline->step(m_cfg.sim_dt);
}

collision::CollisionCoal& PhysicsEngine::collisionManager() {
    if (m_collision_mgr) return *m_collision_mgr;
    throw std::runtime_error("PhysicsEngine::collisionManager(): collision manager not initialized");
}

misc::TimingManager& PhysicsEngine::timings() {
    if (m_timings) return *m_timings;
    throw std::runtime_error("PhysicsEngine::timings(): timings manager not initialized");
}

bool PhysicsEngine::isFinished() const {
    return m_pipeline ? m_pipeline->isFinished() : true;
}

World& PhysicsEngine::world() {
    if (m_world) return *m_world;
    throw std::runtime_error("PhysicsEngine::world(): world not initialized");
}

const World& PhysicsEngine::world() const {
    if (m_world) return *m_world;
    throw std::runtime_error("PhysicsEngine::world() const: world not initialized");
}

void PhysicsEngine::applyForceAt(entt::entity e, const Vector3r& f, const Vector3r& r_local, const Vector3r& tau) {
    if (!m_world) throw std::runtime_error("PhysicsEngine::applyForce(): world not initialized");
    
    const RigidState& state = RigidBody::getState(m_world->ecs(), e);

    Vector3r r_world = state.rotation * r_local; 
    Vector3r induced_tau = r_world.cross(f);
    Vector3r tauEff = tau + induced_tau; 
    Vector3r tauEff_body = state.rotation.transpose() * tauEff;

    m_world->applyForce(e, f, tauEff_body);
}

void PhysicsEngine::disableCollisionBetween(entt::entity a, entt::entity b) {
    if (m_collision_mgr) m_collision_mgr->disablePair(a, b);
}

void PhysicsEngine::setConstraintScalarVelocity(size_t constraintIndex, real_t v) {
    if (constraintIndex < 0 || constraintIndex >= m_world->constraintPatterns().size()) {
        std::cerr << "Warning: setConstraintScalarVelocity: constraintIndex " << constraintIndex << " out of bounds" << std::endl;
        return;
    }
    m_world->constraintPatterns()[constraintIndex]->setScalarVelocity(v);
}

void PhysicsEngine::setConstraintLinearVelocity(size_t constraintIndex, const Vector3r& v) {
    if (constraintIndex < 0 || constraintIndex >= m_world->constraintPatterns().size()) {
        std::cerr << "Warning: setConstraintLinearVelocity: constraintIndex " << constraintIndex << " out of bounds" << std::endl;
        return;
    }
    m_world->constraintPatterns()[constraintIndex]->setLinearVelocity(v);
}

void PhysicsEngine::setConstraintAngularVelocity(size_t constraintIndex, const Vector3r& w) {
    if (constraintIndex < 0 || constraintIndex >= m_world->constraintPatterns().size()) {
        std::cerr << "Warning: setConstraintAngularVelocity: constraintIndex " << constraintIndex << " out of bounds" << std::endl;
        return;
    }
    m_world->constraintPatterns()[constraintIndex]->setAngularVelocity(w);
}
}  // namespace physics
}  // namespace cardillo
