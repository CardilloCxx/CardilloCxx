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

namespace cardillo::solver {

class ConjugateGradientSolver : public SolverBase {
   public:
    explicit ConjugateGradientSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(dyn, cfg) {}

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ConjugateGradient"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    config::Config m_cfg;
    physics::assembly::PgsAssembler m_assembler;
};

}  // namespace cardillo::solver
