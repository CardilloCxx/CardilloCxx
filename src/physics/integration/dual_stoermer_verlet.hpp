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

class DualStoermerVerletSolver : public IntegrationBase {
public:
    explicit DualStoermerVerletSolver(cardillo::World& world, cardillo::physics::DynamicsAssembler& dyn, cardillo::solver::WarmstartProvider* warmstart)
        : IntegrationBase(world, dyn, warmstart) {
        m_dyn.refreshState();
    }

	void step(real_t dt) override;

private:
	// Constraint-space Lagrange multipliers for generalized springs (size = C_dyn)
	VectorXr m_Lambda_g;
};

} // namespace cardillo::integration

