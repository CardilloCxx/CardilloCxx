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

    // Create and own low-level subsystems
    m_timings = std::make_unique<cardillo::misc::TimingManager>();
    m_collision_mgr = std::make_unique<cardillo::collision::CollisionCoal>(*m_world, m_timings.get(), m_cfg);

    // Create pipeline and pass owned subsystem pointers explicitly
    m_pipeline = std::make_unique<cardillo::physics::pipeline::PhysicsPipeline>(*m_world, m_cfg,
                                                                               m_collision_mgr.get(),
                                                                               m_timings.get());
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

cardillo::misc::TimingManager& PhysicsEngine::timings() {
    if (m_timings) return *m_timings;
    throw std::runtime_error("PhysicsEngine::timings(): timings manager not initialized");
}

bool PhysicsEngine::isFinished() const {
    return m_pipeline ? m_pipeline->isFinished() : true;
}

void PhysicsEngine::disableCollisionBetween(entt::entity a, entt::entity b) {
    if (m_collision_mgr) m_collision_mgr->disablePair(a, b);
}

} // namespace physics
} // namespace cardillo
