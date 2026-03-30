#include "physics_pipeline.hpp"

#include <memory>

#include "../assembly/dynamics_assembler.hpp"
#include "../integration/moreau.hpp"
#include "../integration/dual_stoermer_verlet.hpp"
#include "../../io/vtk_writer_binary.hpp"
#include "../../misc/progress/ProgressBar.hpp"

namespace cardillo {
namespace physics {
namespace pipeline {

PhysicsPipeline::~PhysicsPipeline() = default;

PhysicsPipeline::PhysicsPipeline(World& world,
                                 config::Config& cfg,
                                 cardillo::collision::CollisionCoal* collision_mgr,
                                 cardillo::misc::TimingManager* timings,
                                 cardillo::solver::WarmstartProvider* warmstart_provider)
    : m_world(world), m_cfg(cfg), m_collision_mgr(collision_mgr), m_timings(timings), m_warmstart_provider(warmstart_provider) {
    // Create owned components
    m_dyn = std::make_unique<cardillo::physics::DynamicsAssembler>(m_world, m_collision_mgr, m_timings, m_cfg);

    // Choose solver based on config
    using namespace cardillo::solver;
    if (cfg.solver == cardillo::config::SolverType::ProjectedJacobi) {
        m_solver = std::make_unique<ProjectedJacobiSolver>(*m_dyn, m_cfg, m_warmstart_provider);
    } else {
        std::cerr << "Warning: Unrecognized solver type in config; defaulting to ProjectedJacobi.\n";
        m_solver = std::make_unique<ProjectedJacobiSolver>(*m_dyn, m_cfg, m_warmstart_provider);
    }

    // Choose integrator based on config
    using namespace cardillo::integration;
    if (cfg.integrator == cardillo::config::IntegratorType::Moreau) {
        m_integrator = std::make_unique<MoreauSolver>(m_world, *m_solver, *m_dyn, *m_timings, m_cfg);
    } else if (cfg.integrator == cardillo::config::IntegratorType::StoermerVerlet) {
        m_integrator = std::make_unique<DualStoermerVerletIntegrator>(m_world, *m_solver, *m_dyn, *m_timings, m_cfg);
    } else {
        std::cerr << "Warning: Unrecognized integrator type in config; defaulting to Moreau." << std::endl;
        m_integrator = std::make_unique<MoreauSolver>(m_world, *m_solver, *m_dyn, *m_timings, m_cfg);
    }

    // Create VTK writer if output is enabled
    if (cfg.output_interval_steps > 0) {
        m_vtk_writer = std::make_unique<cardillo::io::VtkWriterBinary>(cfg);
    }

    // Create progress bar if sim parameters are available
    if (cfg.sim_dt > 0 && cfg.sim_T > 0) {
        m_total_steps = static_cast<int>(std::ceil(cfg.sim_T / cfg.sim_dt));
        if (m_total_steps > 0) {
            m_pbar = std::make_unique<cardillo::misc::ProgressBar>(static_cast<std::size_t>(m_total_steps), std::cout);
            m_pbar->set_description("Simulating");
        }
    }
}

void PhysicsPipeline::step(real_t dt) {
    if (m_finished) return;

    // Refresh assembler and run integrator
    if (m_dyn) m_dyn->refreshState();
    if (m_integrator) m_integrator->step(dt);

    // Advance internal time/step
    ++m_current_step;
    m_current_time += dt;

    // Output VTK if present
    if (m_vtk_writer) {
        m_vtk_writer->maybeWrite(m_current_step, m_current_time, m_world,
                     m_collision_mgr, m_timings,
                     m_dyn.get());
    }

    // Update progress bar if present
    if (m_pbar) {
        int jorIters = m_solver->lastIterations();
        if (jorIters >= 0) m_pbar->set_postfix("jor=" + std::to_string(jorIters) + "        ");
        m_pbar->update(1);
    }

    if (m_total_steps > 0 && m_current_step >= m_total_steps) {
        m_finished = true;
        if (m_pbar) m_pbar->close();
    }
}

} // namespace pipeline
} // namespace physics
} // namespace cardillo
