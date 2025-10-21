#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"
#include "../solver/projected_jacobi.hpp"
#include "../config/config.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

class MoreauSolver {
public:
	explicit MoreauSolver(cardillo::PhysicsSystem& sys, cardillo::config::Config cfg = cardillo::config::Config{})
			: m_sys(sys), m_dyn(sys), m_pj(m_dyn), m_cfg(std::move(cfg)) {
		m_dyn.assignDofs();
		m_dyn.loadStateFromSystem();
		m_pj.enableWarmStart(true);
		m_pj.setAlpha(m_cfg.pj_alpha);
		m_pj.setCompliance(m_cfg.pj_compliance);
		m_pj.setRelaxation(m_cfg.pj_relaxation);
		m_pj.setMaxIterations(m_cfg.pj_max_iterations);
	}

	// Midpoint rule for unconstrained translation-only point masses
	void stepMidpoint(real_t dt);

	// Access dynamics assembler if needed (e.g., for contacts/W/G)
	cardillo::physics::DynamicsAssembler& dynamics() { return m_dyn; }

private:
	cardillo::PhysicsSystem& m_sys;
	cardillo::physics::DynamicsAssembler m_dyn;
	cardillo::solver::ProjectedJacobiSolver m_pj;
	cardillo::config::Config m_cfg;
};

// No legacy free function; use MoreauSolver directly

}

