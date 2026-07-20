#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <optional>
#include "../../config/config.hpp"
#include "../../misc/types.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/pgs_assembler.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

/// Iterative projected Gauss-Seidel contact solver.
class ProjectedGaussSeidelSolver : public SolverBase {
   public:
    /// Construct the solver from the dynamics assembler and config.
    explicit ProjectedGaussSeidelSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(dyn, cfg) {}

    /// Solve for the corrected velocity increment for one time step.
    VectorXr solve(real_t dt, real_t theta) override;

    /// Solver name used in diagnostics and logs.
    const char* name() const override { return "ProjectedGaussSeidel"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    config::Config m_cfg;
    physics::assembly::PgsAssembler m_assembler;
};

}  // namespace cardillo::solver
