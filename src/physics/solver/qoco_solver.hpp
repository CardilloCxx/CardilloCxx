#pragma once

#include "solver_base.hpp"
#include "../../config/config.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/qoco_assembler.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

class QocoSolver : public SolverBase {
public:
    explicit QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                        const cardillo::config::Config& cfg);
    ~QocoSolver() override;

    VectorXr solve(real_t dt, real_t theta)  override;

    const char* name() const override { return "QocoSolver"; }

private:
    void initQocoSolver(real_t dt, real_t theta, bool first_init = true);
    void updateQocoSolver(real_t dt, real_t theta);

    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::physics::assembly::QocoAssembler m_assembler{m_dyn};
    cardillo::config::Config m_cfg;

    QOCOSolver* m_qoco_solver{nullptr};
};

} // namespace cardillo::solver
