#pragma once

#include <memory>
#include <vector>

#include "../../../../config/config.hpp"
#include "conicxx_assembler.hpp"
#include "../../../assembly/dynamics_assembler.hpp"
#include "../../solver_base.hpp"
#include "../../warmstart.hpp"

#include <conicxx/solver.h>

namespace cardillo::solver {

// Unlike ClarabelSolver/QocoSolver (which rebuild their solver object every step), ConicxxSolver
// keeps one conicxx::Solver for the whole simulation: conicxx::Solver::updateData() when the
// sparsity pattern and the cone structure are unchanged (ConicXX then warm-starts from its previous
// solution), conicxx::Solver::setup() otherwise. Spring and damper rows are compliant zero-cone
// rows by default (conicxx.compliant_rows), see ConicxxAssembler.
class ConicxxSolver : public SolverBase {
   public:
    explicit ConicxxSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~ConicxxSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "ConicxxSolver"; }

   private:
    conicxx::Settings makeSettings() const;
    bool coneSpecChanged(const conicxx::ConeSpec& spec) const;
    // Warm start after a setup() (which discards ConicXX's own): x = current velocities (and
    // percussions), z = previous percussions and matched contact impulses, s = b - Ax + C z_E.
    void setWarmStartFromState(const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const VectorXr& compliance);

    physics::DynamicsAssembler& m_dyn;
    config::Config m_cfg;
    physics::assembly::ConicxxAssembler m_assembler;

    conicxx::Solver m_solver;
    bool m_setup_done{false};
    conicxx::ConeSpec m_last_cone_spec;
};

}  // namespace cardillo::solver
