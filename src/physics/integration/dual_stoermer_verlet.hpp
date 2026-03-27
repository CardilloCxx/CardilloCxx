#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "misc/types.hpp"
#include "integration_base.hpp"
#include "physics/world.hpp"
#include "config/config.hpp"
#include "physics/assembly/dynamics_assembler.hpp"

namespace cardillo::integration {

class DualStoermerVerletIntegrator : public IntegrationBase {
public:
    explicit DualStoermerVerletIntegrator(cardillo::World& world, 
                             cardillo::solver::SolverBase& solver,
                             cardillo::physics::DynamicsAssembler& dyn,
                             cardillo::misc::TimingManager& timings,
                             cardillo::config::Config& config)  : IntegrationBase(world, solver, dyn, timings, config) {}


	void step(real_t dt) override;

private:
	// Constraint-space Lagrange multipliers for generalized springs (size = C_dyn)
	VectorXr m_Lambda_g;
};

} // namespace cardillo::integration

