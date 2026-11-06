#pragma once

#include "../../config/config.hpp"
#include "../../misc/types.hpp"
#include "../assembly/condensed_assembler.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "solver_base.hpp"

namespace cardillo::solver {

// Block-sparse, matrix-free projected Gauss-Seidel on the Delassus blocks of all rows (same rhs()/
// compliance formulas and sign convention as PGS, see condensed_assembler.hpp), built on a
// std::vector<RowBlock>. Sweep modes: sequential Gauss-Seidel or graph-colored (OpenMP-parallel
// within a color); local solve: projection or semismooth Newton; condensed.true_schur eliminates
// the bilateral rows exactly with a block-sparse LDLT in every sweep.
class CondensedSolver : public SolverBase {
   public:
    explicit CondensedSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(dyn, cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;
    const char* name() const override { return "Condensed"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    config::Config m_cfg;
    physics::assembly::CondensedAssembler m_assembler;
};

}  // namespace cardillo::solver
