#pragma once

#include <algorithm>
#include <fstream>

#include "../../../config/config.hpp"
#include "../../assembly/clarabel_assembler.hpp"
#include "../../assembly/dynamics_assembler.hpp"
#include "../solver_base.hpp"
#include "scs_backend.hpp"

namespace cardillo::solver {

/**
 * @brief ADMM contact solver backed by SCS (github.com/cvxgrp/scs) with its direct (AMD + QDLDL)
 * linear-system backend.
 *
 * Solves exactly the cone program ClarabelSolver solves -- the problem data comes from the same
 * ClarabelAssembler (P_mu-scaled contact rows, friction shift in b) -- so results are directly
 * comparable. A fresh scs_init (incl. KKT factorization) is done every step, since A changes with the
 * positions. Unlike the interior-point backends, SCS can be warm started: x from the start-of-step
 * velocities and spring/damper multipliers, y from the contact impulses tracked across steps (the
 * same per-contact storage PGS/PJ warm start from), and optionally the previous step's adapted scale.
 *
 * Only built with -DCARDILLO_WITH_SCS=ON.
 */
class ScsSolver : public SolverBase {
   public:
    explicit ScsSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~ScsSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    real_t lastError() const override { return std::max(m_solver.info().res_pri, m_solver.info().res_dual); }
    const char* name() const override { return "ScsSolver"; }

   private:
    scs::Settings makeSettings() const;
    void fillWarmStart(const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const VectorXr& Smu);
    void writeStats();

    physics::DynamicsAssembler& m_dyn;
    physics::assembly::ClarabelAssembler m_assembler;
    config::Config m_cfg;

    scs::Solver m_solver;
    scs::ConeDims m_dims;
    real_t m_last_scale{-1};

    long m_steps{0};
    long m_inaccurate_steps{0};
    std::ofstream m_stats;
};

}  // namespace cardillo::solver
