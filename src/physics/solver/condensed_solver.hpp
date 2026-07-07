#pragma once

#include "../../config/config.hpp"
#include "../../misc/types.hpp"
#include "../assembly/condensed_assembler.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "solver_base.hpp"

namespace cardillo::solver {

// Block-sparse, matrix-free contact solver generalizing ProjectedGaussSeidel: same rhs()/
// compliance formulas and sign convention (see condensed_assembler.hpp), but built directly on a
// std::vector<RowBlock> instead of Eigen::SparseMatrix, so it can support multiple sweep
// strategies (Jacobi, Gauss-Seidel, graph-colored parallel, chaotic) and local-solve strategies
// (projection, semismooth Newton) as config choices rather than separate classes.
//
// v1 (this file): only "gauss_seidel" sweep mode + "projection" local solve are implemented,
// single-threaded, to validate the block representation against PGS before adding parallelism.
class CondensedSolver : public SolverBase {
   public:
    explicit CondensedSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(dyn, cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;
    const char* name() const override { return "Condensed"; }

   private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::config::Config m_cfg;
    cardillo::physics::assembly::CondensedAssembler m_assembler;
};

}  // namespace cardillo::solver
