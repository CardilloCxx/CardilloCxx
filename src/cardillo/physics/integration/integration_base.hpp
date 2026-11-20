#pragma once

#include <memory>
#include <optional>
#include "../solver/projected_jacobi.hpp"
#include "../solver/solver_base.hpp"
#include "cardillo/misc/types.hpp"
#include "cardillo/physics/assembly/dynamics_assembler.hpp"
#include "cardillo/physics/world.hpp"

namespace cardillo::integration {

/// Observer of the one-step problems (the ffclib export, io/ffclib/ffclib_exporter.hpp). beforeSolve() is called
/// after the step's problem has been assembled, afterSolve() after the solve (and the final check),
/// while the assembler still holds the problem and before the velocities are written back. Observers
/// only read.
class StepObserver {
   public:
    virtual ~StepObserver() = default;
    virtual void beforeSolve(const physics::DynamicsAssembler& dyn) = 0;
    virtual void afterSolve(physics::DynamicsAssembler& dyn, const VectorXr& u, real_t dt, real_t theta, const solver::SolverBase& solver,
                            const std::optional<solver::StepResidual>& residual, real_t solve_seconds) = 0;
};

class IntegrationBase {
   public:
    explicit IntegrationBase(World& world, solver::SolverBase& solver, physics::DynamicsAssembler& dyn, misc::TimingManager& timings,
                             config::Config& config)
        : m_world(world), m_dyn(dyn), m_timings(timings), m_solver(solver), m_config(config) {}

    virtual ~IntegrationBase() = default;
    virtual void step(real_t dt) = 0;

    void explicitPositionUpdate(World& world, real_t h);
    void linearImplicitPositionUpdate(World& world, real_t h);

    /// Final check of the last step (solver.final_check = natural_map); empty if disabled.
    const std::optional<solver::StepResidual>& lastResidual() const { return m_lastResidual; }
    /// Wall-clock time of the solver in the last step (s).
    real_t lastSolveSeconds() const { return m_lastSolveSeconds; }
    /// Observer of every step's problem (not owned; nullptr: none).
    void setStepObserver(StepObserver* observer) { m_observer = observer; }

   protected:
    World& m_world;
    physics::DynamicsAssembler& m_dyn;
    misc::TimingManager& m_timings;
    solver::SolverBase& m_solver;
    config::Config& m_config;
    std::optional<solver::StepResidual> m_lastResidual;
    real_t m_lastSolveSeconds{0};
    StepObserver* m_observer{nullptr};
};

}  // namespace cardillo::integration
