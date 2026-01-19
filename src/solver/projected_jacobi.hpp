#pragma once

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <optional>
#include <algorithm>
#include "../misc/types.hpp"
#include "../physics/dynamics_assembler.hpp"
#include "warmstart.hpp"
#include "../config/config.hpp"

namespace cardillo::solver {

// Simple Projected Jacobi solver for normal contact impulses (percussions)
// Iteration: p_{k+1} = min(0, (R*G - I) p_k + R*b)
// where R = alpha * D^{-1}, D = diag(G), and b = W * v (normal relative speeds)
// Notes:
// - We treat the projection as clamping to non-positive impulses (compressive only).
// - alpha in (0,2) typically for convergence; default alpha=1.
class ProjectedJacobiSolver {
public:
    explicit ProjectedJacobiSolver(cardillo::physics::DynamicsAssembler& dyn,
                                   const cardillo::config::Config& cfg,
                                   WarmstartProvider* cache)
        : m_dyn(dyn)
        , m_alpha(cfg.pj_alpha)
        
        , m_relax(std::clamp<real_t>(cfg.pj_relaxation, 0, 1))
        , m_maxIterations(std::max(1, cfg.pj_max_iterations))
        , m_epsRel(cfg.pj_tol_rel)
        , m_warmStart(cfg.pj_warmstart)
        , m_debug(cfg.debug_pj)
        , m_useNesterov(cfg.pj_nesterov)
        , m_nest_beta_threshold((double)cfg.pj_nesterov_beta_threshold)
        , m_nest_restart_limit(std::max(0, cfg.pj_nesterov_restart_limit))
        , m_wsProvider(cache)
        , m_convCsvDir(cfg.pj_convergence_csv_dir)
    {}

    real_t alpha() const { return m_alpha; }
    int lastIterations() const { return m_lastIterations; }
    real_t lastError() const { return m_lastError; }
    
    // Concatenated API: accept stacked preliminary velocities and return stacked velocities
    VectorXr solve(VectorXr& rhs, real_t tol = 1e-5);

private:
    cardillo::physics::DynamicsAssembler& m_dyn;
    real_t m_alpha{(real_t)1};
    // No internal row-vector warmstart fallback; only via external cache when enabled
    real_t m_relax{(real_t)1};
    int m_maxIterations{1000000};
    bool m_warmStart{true};
    bool m_debug{false};
    bool m_useNesterov{false};
    double m_nest_beta_threshold{0.995};
    int m_nest_restart_limit{4};
    real_t m_epsRel{(real_t)0};
    WarmstartProvider* m_wsProvider{nullptr};
    std::string m_convCsvDir{};
    int m_lastIterations{0};
    real_t m_lastError{(real_t)0};
};

} // namespace cardillo::solver
