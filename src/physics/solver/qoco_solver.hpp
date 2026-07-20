#pragma once

#include "../../config/config.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "../assembly/qoco_assembler.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

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
    void ensureQocoApiLoaded(bool first_init);

    physics::DynamicsAssembler& m_dyn;
    physics::assembly::QocoAssembler m_assembler{m_dyn};
    config::Config m_cfg;

    QOCOSolver* m_qoco_solver{nullptr};

    void* m_qoco_lib_handle{nullptr};
    std::string m_loaded_backend{"uninitialized"};

    using SetDefaultSettingsFn = void (*)(QOCOSettings*);
    using QocoSetupFn = QOCOInt (*)(QOCOSolver*, QOCOInt, QOCOInt, QOCOInt, QOCOCscMatrix*, QOCOFloat*, QOCOCscMatrix*, QOCOFloat*, QOCOCscMatrix*, QOCOFloat*, QOCOInt, QOCOInt, QOCOInt*,
                                    QOCOSettings*);
    using QocoSolveFn = QOCOInt (*)(QOCOSolver*);
    using QocoCleanupFn = QOCOInt (*)(QOCOSolver*);

    SetDefaultSettingsFn m_set_default_settings{nullptr};
    QocoSetupFn m_qoco_setup{nullptr};
    QocoSolveFn m_qoco_solve{nullptr};
    QocoCleanupFn m_qoco_cleanup{nullptr};
};

}  // namespace cardillo::solver
