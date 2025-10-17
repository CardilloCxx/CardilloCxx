#pragma once

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <optional>
#include <algorithm>
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
        : m_dyn(dyn) {
        }

    void setAlpha(real_t a) { m_alpha = a; }
    real_t alpha() const { return m_alpha; }

    // Optional stabilization knobs
    void setCompliance(real_t eps) { m_compliance = std::max<real_t>(0, eps); } // adds eps to each Dii
    void setRelaxation(real_t r) { m_relax = std::clamp<real_t>(r, 0, 1); }     // p_new = (1-r)*p_old + r*proj(...)
    void setMaxIterations(int iters) { m_maxIterations = std::max(1, iters); }
    void enableWarmStart(bool enabled) { m_warmStart = enabled; }
    
    // Convenience overload: accept per-body velocity blocks and flatten internally
    std::vector<VectorXr> iterateWithPreliminaryVelocity(const std::vector<VectorXr>& v_pre_blocks, real_t tol = 1e-5);

private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    real_t m_alpha{(real_t)1};
    std::optional<VectorXr> m_last_p = std::nullopt; // store last result for potential warm start
    real_t m_compliance{(real_t)0};
    real_t m_relax{(real_t)1};
    int m_maxIterations{1000000};
    bool m_warmStart{true};
};

} // namespace cardillo::solver
