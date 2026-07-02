#pragma once

#include <memory>
#include <vector>

#include "../../config/config.hpp"
#include "../assembly/conicxx_assembler.hpp"
#include "../assembly/dynamics_assembler.hpp"
#include "solver_base.hpp"
#include "warmstart.hpp"

#include <conicxx/solver.h>

namespace cardillo::solver {

// Unlike ClarabelSolver/QocoSolver (which have to fully rebuild their solver
// object every step, see ClarabelSolver::updateSolver / QocoSolver::updateQocoSolver),
// ConicxxSolver keeps a single conicxx::Solver alive across the whole
// simulation and uses conicxx::Solver::updateData() for the common case where
// the active contact set (and therefore the sparsity pattern of A and the
// cone structure) hasn't changed since the last step. When it has changed,
// it falls back to a fresh conicxx::Solver::setup() call, exactly as
// documented in conicxx/solver.h's class comment. ConicXX also warm-starts
// itself internally (Settings::warm_start) across successful solve() calls,
// so no explicit setWarmStart() bookkeeping is needed here.
class ConicxxSolver : public SolverBase {
   public:
    explicit ConicxxSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg);
    ~ConicxxSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ConicxxSolver"; }

   private:
    conicxx::Settings makeSettings() const;
    bool coneSpecChanged(const conicxx::ConeSpec& spec) const;

    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::physics::assembly::ConicxxAssembler m_assembler;
    cardillo::config::Config m_cfg;

    conicxx::Solver m_solver;
    bool m_setup_done{false};
    conicxx::ConeSpec m_last_cone_spec;
};

}  // namespace cardillo::solver
