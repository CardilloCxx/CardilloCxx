#pragma once

#include <Eigen/Dense>
#include "../../misc/types.hpp"

namespace cardillo {
namespace solver {

/**
 * @brief Common interface of the contact/constraint solvers (projected Jacobi, condensed, PGS,
 * naive Schur and the optional conic backends) used by the Moreau-theta integrator (see cardillo::integration::MoreauIntegrator).
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

    virtual const char* name() const { return "SolverBase"; }

   protected:
    int m_last_iters{-1};
};

}  // namespace solver
}  // namespace cardillo
