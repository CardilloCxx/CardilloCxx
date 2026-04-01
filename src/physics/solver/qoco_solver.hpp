#pragma once

#include "solver_base.hpp"
#include "../../config/config.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/qoco_assembler.hpp"

namespace cardillo::solver {

class QocoSolver : public SolverBase {
public:
    explicit QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                        const cardillo::config::Config& cfg);

    VectorXr solve(VectorXr& rhs, real_t tol = (real_t)1e-5) override;

    const char* name() const override { return "QocoSolver"; }

private:
    void initQocoSolver();

    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::physics::assembly::QocoAssembler m_assembler{m_dyn};
    cardillo::config::Config m_cfg;

    QOCOSolver* m_qoco_solver{nullptr};
};

} // namespace cardillo::solver
