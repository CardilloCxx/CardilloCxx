#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <optional>
#include "../../config/config.hpp"
#include "../../misc/types.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/pj_assembler.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

// Simple Projected Jacobi solver for normal contact impulses (percussions)
// Iteration: p_{k+1} = -proj_C(R*(W*u_k + b) - p_k)
// where R = alpha * D^{-1}, D = diag(W*Minv*W^T), and C = R_0^- x D(mu*pN) the negated admissible
// set of the contact percussions (pN >= 0, see misc/contact_projection.hpp).
// Notes:
// - alpha in (0,2) typically for convergence.
class ProjectedJacobiSolver : public SolverBase {
   public:
    explicit ProjectedJacobiSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg)
        : m_dyn(dyn), m_assembler(dyn, cfg), m_cfg(cfg) {}

    // Concatenated API: accept stacked preliminary velocities and return stacked velocities
    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ProjectedJacobi"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    physics::assembly::PjAssembler m_assembler;
    config::Config m_cfg;
};

}  // namespace cardillo::solver
