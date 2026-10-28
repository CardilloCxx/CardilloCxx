#pragma once

#include "../../config/config.hpp"
#include "../../misc/types.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/pj_assembler.hpp"
#include "solver_base.hpp"

namespace cardillo::solver {

// Reference solver: eliminates all bilateral rows (springs, dampers, perfect constraints) exactly by
// forming the dense Schur complement of the saddle-point system onto the contact rows, then iterates
// on the contact percussions only.
//
// With the factorized saddle-point matrix S of the projected Jacobi solver (theory manual,
// eq:saddle_point) and the contact force directions W_c (columns), one solve with n_c right-hand sides
// gives X = S^{-1} [W_c; 0; 0]. The velocities and bilateral percussions are affine in the contact
// percussions, x = x_free + X pi_c, and the contact velocities read
//   xi_c = c + G_cc pi_c,   G_cc = W_c^T X_u = W_c^T M^{-1} (W_c - W_b G_bb^{-1} G_bc),
// with G_bb = W_b^T M^{-1} W_b + C_hat and G_bc = W_b^T M^{-1} W_c,
// with c the contact velocity of the free solution (including restitution). Every local contact
// update therefore sees the bilateral reaction exactly, and its step lengths come from the diagonal
// blocks of G_cc instead of the mass-only Delassus blocks used by the condensed solver.
//
// Cost: one sparse LU of S (as projected Jacobi), n_c triangular solves and a dense n_c x n_c
// matrix; each sweep is O(n_c^2). Meant as a reference for the matrix-free solvers, not for large
// scenes.
//
// Contacts are swept by Gauss-Seidel (block by block, step lengths from the diagonal blocks of G,
// pj.alpha * pj.relaxation), optionally with Nesterov momentum and adaptive restart (pj.nesterov).
// Configuration shared with the condensed solver: condensed.local_solve, condensed.projection_rho_strategy,
// condensed.newton_*, and the pj.* keys for tolerances, step length, warm start and Nesterov.
class NaiveSchurSolver : public SolverBase {
   public:
    explicit NaiveSchurSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_assembler(dyn, cfg), m_cfg(cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;
    real_t lastError() const override { return m_lastError; }
    const char* name() const override { return "NaiveSchur"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    physics::assembly::PjAssembler m_assembler;
    config::Config m_cfg;
    real_t m_lastError{(real_t)0};
};

}  // namespace cardillo::solver
