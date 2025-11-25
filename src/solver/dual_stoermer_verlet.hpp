#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"
#include "../solver/projected_jacobi.hpp"
#include "warmstart.hpp"
#include "../config/config.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

class DualStoermerVerletSolver {
public:
		explicit DualStoermerVerletSolver(cardillo::PhysicsSystem& sys)
			: m_sys(sys), m_dyn(sys),
			  m_pj(m_dyn, sys.config(), sys.config().pj_warmstart ? sys.warmstartProvider() : nullptr) {
		// Build initial DOF layout, offsets, and load state
		m_dyn.refreshState();
	}

	// Midpoint rule for unconstrained translation-only point masses
	void stepMidpoint(real_t dt);

	// Access dynamics assembler if needed (e.g., for contacts/W/G)
	cardillo::physics::DynamicsAssembler& dynamics() { return m_dyn; }

private:
	cardillo::PhysicsSystem& m_sys;
	cardillo::physics::DynamicsAssembler m_dyn;
	// warmstart provider is owned by the PhysicsSystem; no local member here
	cardillo::solver::ProjectedJacobiSolver m_pj;
	// Constraint-space Lagrange multipliers for generalized springs (size = C_dyn)
	VectorXr m_Lambda_g;
};

// No legacy free function; use DualStoermerVerletSolver directly

}

