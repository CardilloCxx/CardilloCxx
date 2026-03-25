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
	explicit MoreauSolver(cardillo::World& world)
		: MoreauSolver(world, world.config().moreau_theta) {}

	explicit MoreauSolver(cardillo::World& world, real_t theta)
		: IntegrationBase(world), m_theta(theta) {
		m_dyn.refreshState();
	}

	// Midpoint rule for unconstrained translation-only point masses
	void step(real_t dt) override;
	int lastProjectedJacobiIterations() const override { return m_pj.lastIterations(); }

private:

	// Constraint-space Lagrange multipliers for generalized springs (size = C_dyn)
	VectorXr m_Lambda_g;
	real_t m_theta;
};

} // namespace cardillo::integration

