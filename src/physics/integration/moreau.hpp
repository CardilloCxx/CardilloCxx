#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../solver/projected_jacobi.hpp"
#include "config/config.hpp"
#include "integration_base.hpp"
#include "misc/types.hpp"
#include "physics/assembly/dynamics_assembler.hpp"
#include "physics/world.hpp"

namespace cardillo::integration {

/// Moreau theta integrator for the current rigid-body/contact pipeline.
class MoreauIntegrator : public IntegrationBase {
   public:
    /// Construct the integrator from the shared world, solver, dynamics assembler, and timings.
    explicit MoreauIntegrator(cardillo::World& world, cardillo::solver::SolverBase& solver, cardillo::physics::DynamicsAssembler& dyn, cardillo::misc::TimingManager& timings,
                              cardillo::config::Config& config)
        : IntegrationBase(world, solver, dyn, timings, config) {}

    /// Advance one Moreau theta step.
    void step(real_t dt) override;
};

}  // namespace cardillo::integration
