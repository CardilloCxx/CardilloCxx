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

class ProjectedGaussSeidelSolver : public SolverBase {
   public:
    explicit ProjectedGaussSeidelSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(dyn, cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ProjectedGaussSeidel"; }

   private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::config::Config m_cfg;
    cardillo::physics::assembly::PgsAssembler m_assembler;
};

}  // namespace cardillo::solver
