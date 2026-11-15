#pragma once

#include "../../../../config/config.hpp"
#include "../../../assembly/dynamics_assembler.hpp"
#include "../../solver_base.hpp"
#include "../cone_program_assembler.hpp"
#include "cuopt_backend.hpp"

namespace cardillo::solver {

/**
 * @brief Contact solver backed by NVIDIA cuOpt's GPU barrier method.
 *
 * Solves exactly the cone program ClarabelSolver and ScsSolver solve (same ConeProgramAssembler,
 * P_mu-scaled contact rows, friction shift in b), handed to cuOpt as its Wolfe dual (see
 * cuopt::Model), because cuOpt returns no duals for second-order cones. The problem is assembled on
 * the CPU, copied to the GPU and solved cold every step.
 *
 * Only built with -DCARDILLO_WITH_CUOPT=ON; needs an NVIDIA GPU at run time (the constructor throws
 * without one).
 */
class CuoptSolver : public SolverBase {
   public:
    explicit CuoptSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg);
    ~CuoptSolver() override;

    VectorXr solve(real_t dt, real_t theta) override;

    const char* name() const override { return "CuoptSolver"; }

   private:
    physics::DynamicsAssembler& m_dyn;
    physics::assembly::ConeProgramAssembler m_assembler;
    config::Config m_cfg;

    cuopt::Solver m_solver;
    cuopt::ConeDims m_dims;

    long m_steps{0};
    long m_inaccurate_steps{0};
};

}  // namespace cardillo::solver
