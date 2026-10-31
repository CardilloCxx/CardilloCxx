#pragma once

#include <memory>

#include "../../config/config.hpp"
#include "../../io/csv_writer.hpp"
#include "../world.hpp"

// Forward declarations to avoid heavy includes in the header
namespace cardillo {
namespace physics {
class DynamicsAssembler;
}
}  // namespace cardillo
namespace cardillo {
namespace integration {
class IntegrationBase;
}
}  // namespace cardillo
namespace cardillo {
namespace solver {
class SolverBase;
class WarmstartProvider;
}  // namespace solver
}  // namespace cardillo
namespace cardillo {
namespace io {
class VtkWriter;
}
}  // namespace cardillo
namespace cardillo {
namespace misc {
class ProgressBar;
}
}  // namespace cardillo
namespace cardillo {
namespace collision {
class CollisionCoal;
}
}  // namespace cardillo

namespace cardillo {
namespace physics {
namespace pipeline {

/// Orchestrates collision, assembly, solver, integration, and output for one step.
class PhysicsPipeline {
   public:
    /// Construct the pipeline around a world and engine-owned subsystems.
    PhysicsPipeline(World& world, config::Config& cfg, collision::CollisionCoal* collision_mgr, misc::TimingManager* timings);

    /// Advance the pipeline by one step given dt.
    void step(real_t dt);

    /// Query whether the configured simulation horizon has been reached.
    bool isFinished() const { return m_finished; }

    /// Access the dynamics assembler.
    physics::DynamicsAssembler& dynamicsAssembler() const { return *m_dyn; }
    /// Access the active integrator.
    integration::IntegrationBase& integrator() const { return *m_integrator; }
    /// Access the VTK writer, if output is enabled.
    io::VtkWriter* vtkWriter() const { return m_vtk_writer.get(); }
    /// Access the collision manager.
    collision::CollisionCoal& collisionManager() const { return *m_collision_mgr; }

    ~PhysicsPipeline();

   private:
    World& m_world;
    config::Config& m_cfg;
    collision::CollisionCoal* m_collision_mgr{nullptr};
    misc::TimingManager* m_timings;

    // Owned pipeline components
    std::unique_ptr<physics::DynamicsAssembler> m_dyn;
    std::unique_ptr<integration::IntegrationBase> m_integrator;
    std::unique_ptr<solver::SolverBase> m_solver;
    std::unique_ptr<io::VtkWriter> m_vtk_writer;
    std::unique_ptr<io::CsvWriter> m_tracked_writer;
    // Owned progress bar
    std::unique_ptr<misc::ProgressBar> m_pbar;

    // Internal step/time bookkeeping
    int m_current_step{0};
    real_t m_current_time{0.0};
    int m_total_steps{0};
    bool m_finished{false};
};

}  // namespace pipeline
}  // namespace physics
}  // namespace cardillo
