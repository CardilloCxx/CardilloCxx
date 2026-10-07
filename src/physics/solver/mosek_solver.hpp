#pragma once

#include "../../config/config.hpp"
#include "../assembly/clarabel_assembler.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "mosek_backend.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

/**
 * @brief Interior-point contact solver backed by MOSEK's conic optimizer (Optimizer API).
 *
 * Solves exactly the cone program ClarabelSolver solves -- the problem data comes from the same
 * ClarabelAssembler -- so results are directly comparable. The MOSEK environment (and its license
 * checkout) and task live for the whole simulation: when the dimensions, cone composition and
 * sparsity patterns are unchanged from the previous step only the numerical values are pushed into
 * the existing task, otherwise the task is rebuilt (see mosek::ConicSolver).
 *
 * Only built with -DCARDILLO_WITH_MOSEK=ON.
 */
class MosekSolver : public SolverBase {
   public:
    explicit MosekSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~MosekSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "MosekSolver"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    physics::assembly::ClarabelAssembler m_assembler;
    config::Config m_cfg;

    mosek::ConicSolver m_solver;
    mosek::ConeDims m_dims;
    bool m_first_setup{true};
};

}  // namespace cardillo::solver
