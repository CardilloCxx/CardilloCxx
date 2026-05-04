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

class ProjectedNewtonSolver : public SolverBase {
   public:
    explicit ProjectedNewtonSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg)
        : m_dyn(dyn),
          m_cfg(cfg),
          m_assembler(dyn, cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ProjectedNewton"; }

   private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::physics::assembly::PjAssembler m_assembler;
    cardillo::config::Config m_cfg;
};

}  // namespace cardillo::solver
