#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"
#include "../solver/projected_jacobi.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

class MoreauSolver {
public:
	explicit MoreauSolver(cardillo::PhysicsSystem& sys)
			: m_sys(sys), m_dyn(sys), m_pj(m_dyn) {
		m_dyn.assignDofs();
		m_dyn.loadStateFromSystem();
			m_pj.setAlpha((real_t)0.3);
			m_pj.setCompliance((real_t)1e-6); // tiny compliance to help stability
			m_pj.setRelaxation((real_t)0.9);  // under-relax p updates
			m_pj.setMaxIterations(200000);
			m_pj.enableWarmStart(true);
	}

	// Midpoint rule for unconstrained translation-only point masses
	void stepMidpoint(real_t dt);

	// Access dynamics assembler if needed (e.g., for contacts/W/G)
	cardillo::physics::DynamicsAssembler& dynamics() { return m_dyn; }

private:
	cardillo::PhysicsSystem& m_sys;
	cardillo::physics::DynamicsAssembler m_dyn;
	cardillo::solver::ProjectedJacobiSolver m_pj;
};

// No legacy free function; use MoreauSolver directly

}

