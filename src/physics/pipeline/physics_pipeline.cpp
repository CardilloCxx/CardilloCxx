#include "physics_pipeline.hpp"

#include <memory>

#include "../../collision/collision_coal.hpp"

#include "../../io/csv_writer.hpp"
#include "../../io/vtk_writer_binary.hpp"
#include "../../misc/progress/ProgressBar.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../integration/moreau.hpp"
#include "../solver/clarabel_solver.hpp"
#include "../solver/conicxx_solver.hpp"
#include "../solver/conjugate_gradient.hpp"
#include "../solver/projected_gauss_seidel.hpp"
#include "../solver/projected_jacobi.hpp"
#include "../solver/qoco_solver.hpp"

namespace cardillo {
namespace physics {
namespace pipeline {

PhysicsPipeline::~PhysicsPipeline() = default;

PhysicsPipeline::PhysicsPipeline(cardillo::World& world, cardillo::config::Config& cfg, cardillo::collision::CollisionCoal* collision_mgr, cardillo::misc::TimingManager* timings)
    : m_world(world), m_cfg(cfg), m_collision_mgr(collision_mgr), m_timings(timings) {
    // Use engine-owned collision manager (passed-in)
    // Create owned components
    m_dyn = std::make_unique<cardillo::physics::DynamicsAssembler>(m_world, m_collision_mgr, m_timings, m_cfg);

    // Choose solver based on config
    using namespace cardillo::solver;
    if (cfg.solver == cardillo::config::SolverType::ProjectedJacobi) {
        m_solver = std::make_unique<ProjectedJacobiSolver>(*m_dyn, m_cfg);
    } else if (cfg.solver == cardillo::config::SolverType::ConjugateGradient) {
        m_solver = std::make_unique<ConjugateGradientSolver>(*m_dyn, m_cfg);
    } else if (cfg.solver == cardillo::config::SolverType::ProjectedGaussSeidel) {
        m_solver = std::make_unique<ProjectedGaussSeidelSolver>(*m_dyn, m_cfg);
    } else if (cfg.solver == cardillo::config::SolverType::Qoco) {
        m_solver = std::make_unique<QocoSolver>(*m_dyn, m_cfg);
    } else if (cfg.solver == cardillo::config::SolverType::Clarabel) {
        m_solver = std::make_unique<ClarabelSolver>(*m_dyn, m_cfg);
    } else if (cfg.solver == cardillo::config::SolverType::Conicxx) {
        m_solver = std::make_unique<ConicxxSolver>(*m_dyn, m_cfg);
    } else {
        std::cerr << "Warning: Unrecognized solver type in config; defaulting to ProjectedJacobi.\n";
        m_solver = std::make_unique<ProjectedJacobiSolver>(*m_dyn, m_cfg);
    }

    // Choose integrator based on config
    using namespace cardillo::integration;
    if (cfg.integrator == cardillo::config::IntegratorType::Moreau) {
        m_integrator = std::make_unique<MoreauIntegrator>(m_world, *m_solver, *m_dyn, *m_timings, m_cfg);
    } else {
        std::cerr << "Warning: Unrecognized integrator type in config; defaulting to Moreau." << std::endl;
        m_integrator = std::make_unique<MoreauIntegrator>(m_world, *m_solver, *m_dyn, *m_timings, m_cfg);
    }

    // Create VTK writer if output is enabled
    if (cfg.output_interval_steps > 0) {
        m_vtk_writer = std::make_unique<cardillo::io::VtkWriterBinary>(cfg);
    }

    // Create tracked CSV writer object but do not open file yet; CsvWriter will open lazily on
    // first write
    m_tracked_writer = std::make_unique<cardillo::io::CsvWriter>(m_cfg);

    // Create progress bar if sim parameters are available
    if (cfg.sim_dt > 0 && cfg.sim_T > 0) {
        m_total_steps = static_cast<int>(std::ceil(cfg.sim_T / cfg.sim_dt));
        if (m_total_steps > 0) {
            m_pbar = std::make_unique<cardillo::misc::ProgressBar>(static_cast<std::size_t>(m_total_steps), std::cout);
            m_pbar->set_description("Simulating");
        }
    }
}

/**
 * @brief Advances the simulation by one time step of size dt. This method is responsible for:
 * 1. Refreshing the dynamics assembler state
 * 2. Running the integrator to advance the simulation state
 * 3. Writing output (VTK, tracked CSV) if enabled
 * 4. Updating the progress bar if enabled
 * 5. Tracking the current simulation time and step, and setting the finished flag when done
 * @param dt The time step size to advance the simulation by
 * @note This method may be called multiple times in a loop to advance the simulation through time
 */
void PhysicsPipeline::step(real_t dt) {
    if (m_finished) return;

    if (m_current_step == 0 && m_vtk_writer) {
        m_vtk_writer->maybeWrite(m_current_step, m_current_time, m_world, m_collision_mgr, m_timings, m_dyn.get());
    }

    // Run integrator (it refreshes the assembler state itself at the start of its step)
    if (m_integrator) m_integrator->step(dt);

    // Advance internal time/step
    ++m_current_step;
    m_current_time += dt;

    // Output VTK if present
    if (m_vtk_writer) {
        m_vtk_writer->maybeWrite(m_current_step, m_current_time, m_world, m_collision_mgr, m_timings, m_dyn.get());
    }

    // Write tracked state to CSV
    if (m_tracked_writer) {
        m_tracked_writer->writeTrackedState(m_current_time, m_world.ecs());
    }

    // Update progress bar if present
    if (m_pbar) {
        int iters = m_solver->lastIterations();
        if (iters >= 0) m_pbar->set_postfix("iters=" + std::to_string(iters) + "        ");
        m_pbar->update(1);
    }

    if (m_total_steps > 0 && m_current_step >= m_total_steps) {
        m_finished = true;
        if (m_pbar) m_pbar->close();
    }
}

}  // namespace pipeline
}  // namespace physics
}  // namespace cardillo
