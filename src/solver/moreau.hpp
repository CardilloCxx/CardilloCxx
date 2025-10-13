#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

class MoreauSolver {
public:
	explicit MoreauSolver(cardillo::PhysicsSystem& sys)
		: m_sys(sys), m_dyn(sys) {
		m_dyn.assignDofs();
		m_dyn.loadStateFromSystem();
	}

	// Midpoint rule for unconstrained translation-only point masses
	void stepMidpoint(real_t dt);

	// Access dynamics assembler if needed (e.g., for contacts/W/G)
	cardillo::physics::DynamicsAssembler& dynamics() { return m_dyn; }

private:
	cardillo::PhysicsSystem& m_sys;
	cardillo::physics::DynamicsAssembler m_dyn;
};

// No legacy free function; use MoreauSolver directly

}

