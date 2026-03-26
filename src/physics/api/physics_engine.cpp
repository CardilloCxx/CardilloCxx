#include "physics_engine.hpp"

#include "../../collision/collision_coal.hpp"
#include "../misc/timings/TimingManager.hpp"
#include "../../solver/warmstart.hpp"
#include "../pipeline/physics_pipeline.hpp"

namespace cardillo {
namespace physics {

PhysicsEngine::PhysicsEngine() : m_world(std::make_unique<World>()) {}

PhysicsEngine::PhysicsEngine(const config::Config& cfg) : PhysicsEngine() {
    initFromConfig(cfg);
}

PhysicsEngine::PhysicsEngine(std::unique_ptr<cardillo::World> world, const config::Config& cfg)
    : m_world(std::move(world)) {
    initFromConfig(cfg);
}

void PhysicsEngine::initFromConfig(const config::Config& cfg) {
    m_cfg = cfg;
    if (!m_world) m_world = std::make_unique<World>();
    m_world->setConfig(m_cfg);

    // Create and own low-level subsystems
    m_collision_mgr = std::make_unique<cardillo::collision::CollisionCoal>();
    m_collision_mgr->registerSystem(m_world.get());
    m_timings = std::make_unique<cardillo::misc::TimingManager>();
    m_warmstart_provider = std::make_unique<cardillo::solver::WarmstartCache>();

    // Inform world about non-owning pointers to these subsystems
    m_world->setCollisionManager(m_collision_mgr.get());
    m_world->setTimings(m_timings.get());
    m_world->setWarmstartProvider(m_warmstart_provider.get());

    // Create pipeline
    m_pipeline = std::make_unique<cardillo::physics::pipeline::PhysicsPipeline>(*m_world, m_cfg,
                                                                               m_collision_mgr.get(),
                                                                               m_timings.get(),
                                                                               m_warmstart_provider.get());
}

void PhysicsEngine::step(real_t dt) {
    if (m_pipeline) m_pipeline->step(dt);
}

void PhysicsEngine::step() {
    real_t dt = (m_cfg.sim_dt > (real_t)0) ? m_cfg.sim_dt : m_world->config().sim_dt;
    if (m_pipeline) m_pipeline->step(dt);
}

collision::CollisionCoal& PhysicsEngine::collisionManager() {
    if (m_collision_mgr) return *m_collision_mgr;
    return m_world->collisionManager();
}

cardillo::misc::TimingManager& PhysicsEngine::timings() {
    if (m_timings) return *m_timings;
    return m_world->timings();
}

bool PhysicsEngine::isFinished() const {
    return m_pipeline ? m_pipeline->isFinished() : true;
}

} // namespace physics
} // namespace cardillo
