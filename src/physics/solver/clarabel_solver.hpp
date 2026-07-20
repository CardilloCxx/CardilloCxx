#pragma once

#include <memory>
#include <vector>

#include "../../config/config.hpp"
#include "../assembly/clarabel_assembler.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

#include <clarabel.hpp>

namespace cardillo::solver {

class ClarabelSolver : public SolverBase {
   public:
    explicit ClarabelSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~ClarabelSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ClarabelSolver"; }

   private:
    void initSolver(real_t dt, real_t theta, bool first_init = true);
    void updateSolver(real_t dt, real_t theta);

    physics::DynamicsAssembler& m_dyn;
    physics::assembly::ClarabelAssembler m_assembler;
    config::Config m_cfg;

    std::unique_ptr<clarabel::DefaultSolver<double>> m_solver;
};

}  // namespace cardillo::solver
