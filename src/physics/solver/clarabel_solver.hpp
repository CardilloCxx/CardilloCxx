#pragma once

#include <memory>
#include <vector>

#include "solver_base.hpp"
#include "../../config/config.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/clarabel_assembler.hpp"
#include "warmstart.hpp"

#include <clarabel.hpp>

namespace cardillo::solver {

class ClarabelSolver : public SolverBase {
public:
    explicit ClarabelSolver(cardillo::physics::DynamicsAssembler& dyn,
                            const cardillo::config::Config& cfg);
    ~ClarabelSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ClarabelSolver"; }

private:
    void initSolver(real_t dt, real_t theta, bool first_init = true);
    void updateSolver(real_t dt, real_t theta);

    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::physics::assembly::ClarabelAssembler m_assembler;
    cardillo::config::Config m_cfg;

    std::unique_ptr<clarabel::DefaultSolver<double>> m_solver;
};

} // namespace cardillo::solver
