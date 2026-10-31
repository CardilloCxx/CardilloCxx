#pragma once

#include "../../../../config/config.hpp"
#include "../../../assembly/dynamics_assembler.hpp"
#include "qoco_assembler.hpp"
#include "../../solver_base.hpp"
#include "../../warmstart.hpp"

namespace cardillo::solver {

class QocoSolver : public SolverBase {
   public:
    explicit QocoSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~QocoSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "QocoSolver"; }

   private:
    void initQocoSolver(real_t dt, real_t theta, bool first_init = true);
    void updateQocoSolver(real_t dt, real_t theta);

    physics::DynamicsAssembler& m_dyn;
    physics::assembly::QocoAssembler m_assembler{m_dyn};
    config::Config m_cfg;

    QOCOSolver* m_qoco_solver{nullptr};
};

}  // namespace cardillo::solver
