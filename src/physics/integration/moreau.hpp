#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "misc/types.hpp"
#include "integration_base.hpp"
#include "physics/world.hpp"
#include "solver/projected_jacobi.hpp"
#include "solver/warmstart.hpp"
#include "config/config.hpp"
#include "physics/assembly/dynamics_assembler.hpp"

namespace cardillo::integration {

class MoreauSolver : public IntegrationBase {
public:
    explicit MoreauSolver(cardillo::World& world, 
                             cardillo::solver::SolverBase& solver,
                             cardillo::physics::DynamicsAssembler& dyn,
                             cardillo::misc::TimingManager& timings,
                             cardillo::config::Config& config)  : IntegrationBase(world, solver, dyn, timings, config) {}

	// Midpoint rule for unconstrained translation-only point masses
	void step(real_t dt) override;

private:

	// Constraint-space Lagrange multipliers for generalized springs (size = C_dyn)
	VectorXr m_Lambda_g;
	real_t m_theta;
};

} // namespace cardillo::integration

