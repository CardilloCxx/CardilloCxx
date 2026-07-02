#pragma once

#include <Eigen/Dense>
#include "../../misc/types.hpp"

namespace cardillo {
namespace solver {

/**
 * @brief Common interface for the velocity-level contact/constraint solvers (Projected
 * Jacobi/Gauss-Seidel, Conjugate Gradient, and the QOCO/Clarabel cone-program solvers) used by
 * the Moreau-Jean theta-method integrator (see cardillo::integration::MoreauIntegrator).
 */
class SolverBase {
   public:
    virtual ~SolverBase() = default;

    /**
     * @brief Solves one time step's velocity-level system and returns the corrected body
     * velocities (stacked per-body, matching DynamicsAssembler::bodyVelOffsets()).
     * @param dt Time step size; implementations assume `dt > 0` and do not validate it.
     * @param theta Moreau-Jean theta-method splitting parameter; implementations assume
     * `theta` is in `(0, 1]` and do not validate it. A zero/invalid value will silently produce
     * inf/nan through downstream divisions rather than raising an error.
     */
    virtual VectorXr solve(real_t dt, real_t theta) = 0;

    // Diagnostics
    virtual int lastIterations() const { return m_last_iters; }
    virtual real_t lastError() const { return (real_t)0; }

    virtual const char* name() const { return "SolverBase"; }

   protected:
    int m_last_iters{-1};
};

}  // namespace solver
}  // namespace cardillo
