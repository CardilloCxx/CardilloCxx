#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../solver/projected_jacobi.hpp"
#include "cardillo/config/config.hpp"
#include "integration_base.hpp"
#include "cardillo/misc/types.hpp"
#include "cardillo/physics/assembly/dynamics_assembler.hpp"
#include "cardillo/physics/world.hpp"

namespace cardillo::integration {

/// Moreau theta integrator for the current rigid-body/contact pipeline.
class MoreauIntegrator : public IntegrationBase {
   public:
    /// Construct the integrator from the shared world, solver, dynamics assembler, and timings.
    explicit MoreauIntegrator(World& world, solver::SolverBase& solver, physics::DynamicsAssembler& dyn, misc::TimingManager& timings,
                              config::Config& config)
        : IntegrationBase(world, solver, dyn, timings, config) {}

    /// Advance one Moreau theta step.
    void step(real_t dt) override;
};

}  // namespace cardillo::integration
