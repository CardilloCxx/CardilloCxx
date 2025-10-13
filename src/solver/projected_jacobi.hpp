#pragma once

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <optional>
#include "../misc/types.hpp"
#include "../physics/dynamics_assembler.hpp"

namespace cardillo::solver {

// Simple Projected Jacobi solver for normal contact impulses (percussions)
// Iteration: p_{k+1} = min(0, (R*G - I) p_k + R*b)
// where R = alpha * D^{-1}, D = diag(G), and b = W * v (normal relative speeds)
// Notes:
// - We treat the projection as clamping to non-positive impulses (compressive only).
// - alpha in (0,2) typically for convergence; default alpha=1.
class ProjectedJacobiSolver {
public:
    explicit ProjectedJacobiSolver(cardillo::physics::DynamicsAssembler& dyn)
        : m_dyn(dyn) {}

    void setAlpha(real_t a) { m_alpha = a; }
    real_t alpha() const { return m_alpha; }

    // Perform N iterations starting from optional initial p0 (else zeros)
    // Returns the resulting impulse vector p
    VectorXr iterate(std::optional<RefVectorXr> p0 = std::nullopt, real_t tol = 1e-5);

    // Iteration using a provided preliminary velocity v_pre for b = W * v_pre
    VectorXr iterateWithPreliminaryVelocity(const VectorXr& v_pre, std::optional<RefVectorXr> p0 = std::nullopt, real_t tol = 1e-5);

private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    real_t m_alpha{(real_t)1};
};

} // namespace cardillo::solver
