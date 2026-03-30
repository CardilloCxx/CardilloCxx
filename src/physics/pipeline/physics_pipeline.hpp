#pragma once

#include <memory>

#include "../world.hpp"
#include "../../io/csv_writer.hpp"
#include "../../config/config.hpp"

// Forward declarations to avoid heavy includes in the header
namespace cardillo { namespace physics { class DynamicsAssembler; } }
namespace cardillo { namespace integration { class IntegrationBase; } }
namespace cardillo { namespace solver { class SolverBase; class WarmstartProvider; } }
namespace cardillo { namespace io { class VtkWriterBinary; } }
namespace cardillo { namespace misc { class ProgressBar; } }

namespace cardillo {
namespace physics {
namespace pipeline {

class PhysicsPipeline {
public:
    // Construct pipeline with references/pointers to the world and owned subsystems
    PhysicsPipeline(World& world,
                    config::Config& cfg,
                    cardillo::collision::CollisionCoal* collision_mgr,
                    cardillo::misc::TimingManager* timings);

    // Advance the pipeline by one step given dt; pipeline tracks current step/time internally
    void step(real_t dt);

    // Query finished state
    bool isFinished() const { return m_finished; }

    // Getters for owned components
    cardillo::physics::DynamicsAssembler& dynamicsAssembler() { return *m_dyn; }
    cardillo::integration::IntegrationBase& integrator() { return *m_integrator; }
    cardillo::io::VtkWriterBinary* vtkWriter() { return m_vtk_writer.get(); }
    

    ~PhysicsPipeline();

private:
    World& m_world;
    config::Config& m_cfg;
    cardillo::collision::CollisionCoal* m_collision_mgr;
    cardillo::misc::TimingManager* m_timings;
    std::unique_ptr<cardillo::solver::WarmstartProvider> m_warmstart_provider;
    
    
    // Owned pipeline components
    std::unique_ptr<cardillo::physics::DynamicsAssembler> m_dyn;
    std::unique_ptr<cardillo::integration::IntegrationBase> m_integrator;
    std::unique_ptr<cardillo::solver::SolverBase> m_solver;
    std::unique_ptr<cardillo::io::VtkWriterBinary> m_vtk_writer;
    std::unique_ptr<cardillo::io::CsvWriter> m_tracked_writer;
    // Owned progress bar
    std::unique_ptr<cardillo::misc::ProgressBar> m_pbar;

    // Internal step/time bookkeeping
    int m_current_step{0};
    real_t m_current_time{0.0};
    int m_total_steps{0};
    bool m_finished{false};
};

} // namespace pipeline
} // namespace physics
} // namespace cardillo
