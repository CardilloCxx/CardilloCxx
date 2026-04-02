#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "misc/types.hpp"
#include "integration_base.hpp"
#include "physics/world.hpp"
#include "../solver/projected_jacobi.hpp"
#include "config/config.hpp"
#include "physics/assembly/dynamics_assembler.hpp"

namespace cardillo::integration {

class MoreauIntegrator : public IntegrationBase {
public:
    explicit MoreauIntegrator(cardillo::World& world, 
                             cardillo::solver::SolverBase& solver,
                             cardillo::physics::DynamicsAssembler& dyn,
                             cardillo::misc::TimingManager& timings,
                             cardillo::config::Config& config)  : IntegrationBase(world, solver, dyn, timings, config) {}

	// Midpoint rule for unconstrained translation-only point masses
	void step(real_t dt) override;
};

} // namespace cardillo::integration

