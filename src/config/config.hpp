#include <optional>
#pragma once

#include <string>
#include "misc/types.hpp"

namespace cardillo::config {

enum class IntegratorType { Moreau };
enum class SolverType { ProjectedJacobi, ConjugateGradient, ProjectedGaussSeidel, Qoco, Clarabel, Conicxx, Condensed };

struct Config {
    // Simulation settings
    real_t sim_T{(real_t)5.0};    // sim.T - total simulation time
    real_t sim_dt{(real_t)1e-3};  // sim.dt - time step
    // Gravity vector (world frame) used by examples and World
    Vector3r sim_gravity{(real_t)0, (real_t)0, (real_t)-9.81};  // sim.gravity = gx gy gz

    // Output settings
    int output_interval_steps{50};                // output.interval_steps - VTK write interval in steps
    std::string output_folder{"./vtk_out"};       // output.folder
    std::string output_filename_prefix{"scene"};  // output.filename_prefix
    bool output_write_contact_manifolds{true};    // output.write_contact_manifolds

    // Collision settings
    // collision.broadphase: one of [dynamic_aabb, dynamic_aabb_array, naive, sap, ssap]
    std::string collision_broadphase{"dynamic_aabb"};
    // Disable all collisions globally (no collidable components, no contact computation)
    bool collision_disable_all{false};
    // Security margin for collision request (>=0). Small positive helps robustness.
    real_t collision_security_margin{(real_t)0.0};  // collision.security_margin
    // Maximum raw contacts requested per pair from narrowphase
    int collision_max_raw_contacts{1024};  // collision.max_raw_contacts
    // Maximum number of contact patches requested per pair
    std::size_t collision_max_patches{4};  // collision.max_patches
    // Prefer using patch vertices over raw contacts when available
    bool collision_use_patch_vertices{true};  // collision.use_patch_vertices
    // Maximum distance (in meters) to consider a current contact matched to a previous one
    // Used for warmstarting across steps. Set via `collision.match_max_dist`.
    real_t collision_match_max_dist{(real_t)0.02};
    // Minimum separation between contacts within the same entity-entity pair (meters).
    // Contacts closer than this will be deduplicated (keeping deeper penetrations).
    // Set via `collision.min_pair_contact_distance`. Set to 0 to disable.
    real_t collision_min_pair_contact_distance{(real_t)0.0};
    // Friction settings
    bool friction_enable{false};              // friction.enable - extend W with tangents
    std::string friction_combine{"min"};      // friction.combine - one of [min, arith, geom]
    real_t friction_default_mu{(real_t)0.3};  // friction.default_mu - default C_Friction when not set on entities

    // Restitution settings
    real_t restitution_default_normal{(real_t)0.0};      // restitution.default_normal
    real_t restitution_default_tangential{(real_t)0.0};  // restitution.default_tangential

    // Debug/diagnostics
    bool debug_rb{false};    // debug.rb  - enable rigid-body contact/W diagnostics in Moreau
    bool debug_pj{false};    // debug.pj  - enable ProjectedJacobi iteration logging
    bool debug_mesh{false};  // debug.mesh - print mesh normalization info (volume, COM, inertia)

    // Solver selection
    SolverType solver{SolverType::ProjectedJacobi};
    int pj_max_iterations{200000};                     // pj.max_iterations
    real_t pj_tol_abs{(real_t)1e-5};                   // pj.tol_abs
    real_t pj_tol_rel{(real_t)1e-5};                   // pj.tol_rel
    // WARNING: these two keys mean genuinely different things in PJ vs. PGS/condensed -- not just
    // different tuning, a structurally different formula. PGS/condensed multiply BOTH onto the
    // per-iteration delta every sweep: dlambda = relaxation * GiiInv * r, with an EXTRA *alpha
    // applied only to contact rows (so relaxation is a uniform SOR-style damping factor across
    // every row kind, alpha an additional contact-only one). PJ never reads pj_relaxation at all
    // (grepped: zero references in projected_jacobi.cpp) and instead bakes alpha directly into the
    // contact preconditioner ONCE per solve() call: R_block = alpha * D^{-1} (see
    // buildBlockPreconditioner()/rdiag_sparse() in projected_jacobi.cpp). Mathematically alpha's
    // role there matches Anitescu-Tasora's cone-complementarity framing: a step size bounded by the
    // reciprocal of the contact Delassus operator's Lipschitz constant, not an arbitrary damping
    // knob -- confirmed empirically for condensed too (see CONDENSED_SOLVER_REPORT.md's alpha/
    // relaxation review): the power-iteration spectral-radius estimate already built for
    // pj.chebyshev, run at alpha=1, crosses 1 at almost exactly the alpha domino's own config
    // comments arrived at by hand-tuning for condensed.sweep_mode=jacobi. Neither `alpha` nor
    // `relaxation` is currently derived from that estimate automatically for either solver --
    // hand-tuned per scene today. See the report for the full theory-vs-implementation writeup and
    // the roadmap of what's untested (e.g. relaxation > 1, true over-relaxation, has never been
    // tried anywhere in this codebase's example scenes despite being provably safe for the
    // symmetric-positive-definite bilateral sub-system per the Ostrowski-Reich theorem).
    real_t pj_relaxation{(real_t)1.0};                 // pj.relaxation
    real_t pj_alpha{(real_t)0.3};                      // pj.alpha
    bool pj_nesterov{false};                           // pj.nesterov (enable Nesterov acceleration)
    real_t pj_nesterov_beta_threshold{(real_t)1.0};    // pj.nesterov_beta_threshold
    int pj_nesterov_restart_limit{10000};              // pj.nesterov_restart_limit
    // Chebyshev semi-iterative acceleration: an alternative to Nesterov, estimating the iteration
    // operator's spectral radius via power iteration each solve() call. Mutually exclusive with
    // pj_nesterov (pj_nesterov takes precedence if both are set). Shared by both
    // solver.type=projected_jacobi (matrix-based power iteration through PjAssembler::solveS())
    // and solver.type=condensed (matrix-free power iteration reusing condensed's own
    // block-sweep functions with projection/Newton disabled and a zero rhs -- see
    // estimateSpectralRadius() in condensed_solver.cpp). Verified on three scenes against a true
    // (non-accelerated) baseline: 3.4-10.7x fewer sweeps than no acceleration, matching or beating
    // pj_nesterov in every case tested (up to ~1.8x fewer sweeps than Nesterov on dense
    // rigid-contact scenes; within ~3% of Nesterov, still both a ~3.4x win over baseline, on an
    // elastic+frictional scene) -- for PJ. On condensed, ~20% fewer sweeps than no acceleration and
    // competitive with/slightly ahead of Nesterov on both domino and slinky (see
    // docs/chapters/solvers/condensed.rst). Not yet default -- see the accompanying report for the
    // full methodology before relying on these numbers for a scene not yet tested.
    bool pj_chebyshev{false};  // pj.chebyshev
    // Anderson acceleration (Type-II / unconstrained-least-squares reformulation): mixes the last
    // pj_anderson_m residual differences via a small least-squares solve each sweep. Precedence
    // when multiple acceleration flags are set: pj_nesterov > pj_chebyshev > pj_anderson.
    // Window size matters a lot in practice -- see the accompanying report's benchmark sweep
    // across window sizes before picking one for a new scene; m=1 (which has a closed-form
    // solution, no least-squares solve needed) was found not competitive on its own.
    bool pj_anderson{false};  // pj.anderson
    int pj_anderson_m{2};     // pj.anderson_m (history window size, >= 1; 2 was the best of {1,2,3,5,8,12,20} across all three benchmark scenes -- see the report)
    bool pj_warmstart{true};  // pj.warmstart (enable warmstart & cache)
    bool pj_rdiag_true_delassus{false};                // pj.rdiag_true_delassus
    std::string pj_convergence_csv_dir{""};            // pj.convergence_csv_dir (empty to disable)

    // Condensed (matrix-free block-relaxation) solver settings. Iteration/tolerance/relaxation/
    // warmstart/debug knobs are shared with PJ/PGS -- see pj_max_iterations, pj_tol_abs,
    // pj_tol_rel, pj_relaxation, pj_alpha, pj_warmstart, debug_pj above.
    std::string condensed_sweep_mode{"gauss_seidel"};    // condensed.sweep_mode: jacobi | gauss_seidel | colored | chaotic
    std::string condensed_local_solve{"projection"};     // condensed.local_solve: projection | newton
    std::string condensed_ordering{"natural"};           // condensed.ordering: natural | dominant
    int condensed_num_threads{0};                        // condensed.num_threads (0 = OpenMP default)
    int condensed_newton_max_iters{8};                   // condensed.newton_max_iters
    real_t condensed_newton_tol{(real_t)1e-10};          // condensed.newton_tol
    // How the per-contact Newton step's diagonal preconditioner is derived from the block's own
    // local Delassus matrix -- "split" (default, this codebase's original strategy: normal and
    // tangential treated as separate 1- and 2-dof sub-problems) or "full" (a single rho from the
    // FULL 3x3 block's largest eigenvalue, capturing normal-tangential coupling that "split"
    // ignores). Matches Siconos's own compute_rho_split_spectral_norm/compute_rho_spectral_norm
    // (fc3d_AlartCurnier_functions.c) respectively -- checked directly against that codebase. See
    // CONDENSED_SOLVER_REPORT.md for a head-to-head comparison before switching a scene to "full".
    std::string condensed_newton_rho_strategy{"split"}; // condensed.newton_rho_strategy: split | full
    // How GiiInv (the plain, non-Newton per-block preconditioner used by every sweep mode's
    // fallback/default path) is derived for ContactFrictional (3-row, normal+2-tangential) blocks
    // only -- "diagonal" (default) is the existing per-DOF-independent inverse
    // (1/Gii(i,i) each), deliberately ignoring intra-block coupling (matches PgsAssembler::
    // DinvDiag(); using the FULL block inverse here was tried and found to diverge on
    // strongly-coupled multi-contact systems like domino -- see condensed_assembler.cpp's
    // updateCompliance()). "split"/"full" instead scale the diagonal by the SAME closed-form
    // rho estimate (misc/contact_rho.hpp) already computed for the Newton local solve --
    // strictly more conservative than the full inverse (1/lambda_max <= 1/any diagonal entry for
    // a real symmetric matrix), while still capturing tangential-tangential (and, for "full",
    // normal-tangential) coupling the plain diagonal ignores entirely. Computed once per step in
    // updateCompliance() (Gii is fixed for the whole step's sweep), not per sweep iteration. See
    // CONDENSED_SOLVER_REPORT.md for a head-to-head comparison before switching a scene off
    // "diagonal".
    std::string condensed_projection_rho_strategy{"diagonal"}; // condensed.projection_rho_strategy: diagonal | split | full
    int condensed_chaotic_reshuffle_interval{50};        // condensed.chaotic_reshuffle_interval (sweeps between reshuffles)
    unsigned condensed_chaotic_seed{12345u};             // condensed.chaotic_seed
    // condensed.true_schur: exactly eliminate the bilateral (spring+damper) rows every outer
    // iteration via a block-sparse LDLT factorization (see BlockSparseLDLT, CondensedAssembler::
    // buildBilateralFactorization()), instead of letting them participate in the Gauss-Seidel/
    // Jacobi/colored/chaotic sweep like every other row. Default false -- preserves today's
    // validated behavior; only affects scenes with springs/dampers (a no-op on e.g. domino).
    bool condensed_true_schur{false};                    // condensed.true_schur

    // EXPERIMENTAL, opt-in, default off: derive `alpha` for THIS solve() call from the same
    // zero-rhs/linearOnly power iteration built for pj.chebyshev, instead of using the fixed
    // pj_alpha every scene currently hand-tunes once and leaves static for the whole run. Measures
    // the raw (unclamped) spectral radius at alpha=1, then picks alpha so the resulting operator's
    // estimated radius is approximately condensed_auto_alpha_target_rho (a linear-regime
    // extrapolation: rho(alpha) ~= alpha*(rho(1)+1) - 1 for alpha large enough to be dominated by
    // the largest eigenvalue, confirmed empirically on domino -- see CONDENSED_SOLVER_REPORT.md).
    // Recomputed fresh every step, unlike a static hand-tuned alpha -- motivated by the measured
    // finding that the true spectral radius is state-dependent and can swing between stable and
    // divergent within a few steps of the same run. Not yet validated widely enough to change the
    // default; see the report for head-to-head comparisons against hand-tuned alpha before relying
    // on this for a scene not yet tested there.
    bool condensed_auto_alpha{false};                    // condensed.auto_alpha
    real_t condensed_auto_alpha_target_rho{(real_t)0.7}; // condensed.auto_alpha_target_rho

    // EXPERIMENTAL, opt-in, default off, only meaningful with pj.chebyshev: instead of the
    // separate zero-rhs/linearOnly power-iteration estimate above (measured to be overly
    // pessimistic for contact-rich scenes -- it ignores the friction-cone projection's own
    // stabilizing effect entirely), run a few real (unaccelerated, actually-projected) sweeps
    // first and estimate rho from the OBSERVED residual-decay ratio instead (classical "adaptive
    // SOR/Chebyshev", e.g. Manteuffel 1977/Ashby 1985).
    // MEASURED RESULT (see CONDENSED_SOLVER_REPORT.md): on domino, this performed WORSE than the
    // upfront linear estimate (1847 vs 1415 total sweeps), not better. Tracing the observed ratio
    // per warmup iteration found why: it does not settle to a stable value even over 20
    // iterations -- drifting from ~0.3 up through ~0.9 and transiently EXCEEDING 1 (residual
    // growing) as the active contact set keeps resolving. The classical technique assumes a
    // linearly-convergent iteration with a well-defined asymptotic ratio; this codebase's
    // Gauss-Seidel/Jacobi-with-projection sweep doesn't obviously have one on any practically
    // affordable observation window. Left available (not removed -- it's not wrong, just not
    // beneficial so far) in case a longer/smarter ratio estimator or a different scene changes
    // this finding; do not enable without re-measuring on your own scene first.
    bool condensed_chebyshev_adaptive_rho{false};        // condensed.chebyshev_adaptive_rho

    real_t constraint_bias_factor{(real_t)0.001};  // constraint_bias_factor - Baumgarte-style bias factor for position error correction (0 to disable)

    // Interior-point solver settings shared across QOCO, Clarabel and ConicXX
    std::string qoco_backend{"auto"};         // qoco.backend [auto, cpu, cuda]
    real_t ip_kkt_static_reg{(real_t)1e-8};   // solver.kkt_static_reg
    real_t ip_kkt_dynamic_reg{(real_t)1e-8};  // solver.kkt_dynamic_reg
    int ip_iter_ref_iters{2};                 // solver.iter_ref_iters

    // ConicXX-only: unlike QOCO/Clarabel it can reuse its KKT factorization and
    // warm-start its iterate across steps when the active contact set hasn't
    // changed (see docs/chapters/solvers/interior_point.rst).
    bool conicxx_warm_start{true};  // conicxx.warm_start

    // ConicXX-only: which sparse LDL^T backend factorizes/solves the KKT system.
    // One of "eigen" (Eigen::SimplicialLDLT, always available), "qdldl" (the
    // backend QOCO/Clarabel use -- ConicXX's own default), or
    // "regularized_ldlt" (a custom Davis/ECOS-style backend correcting bad
    // pivots inline during elimination instead of KktSystem's outer
    // refactorize-from-scratch retry loop). "qdldl"/"regularized_ldlt" fall
    // back to "eigen" with a one-time warning if ConicXX wasn't built with
    // CONICXX_WITH_QDLDL.
    std::string conicxx_linear_solver{"qdldl"};  // conicxx.linear_solver [eigen, qdldl, regularized_ldlt]

    // ConicXX-only: remaining conicxx::Settings fields not already covered by
    // the shared ip_kkt_*/pj_tol_*/pj_max_iterations fields above.
    real_t conicxx_dynamic_reg_eps{(real_t)1e-14};      // conicxx.dynamic_reg_eps -- pivot-magnitude threshold triggering a regularization bump
    real_t conicxx_refine_tol{(real_t)1e-12};           // conicxx.refine_tol -- target relative residual for iterative refinement
    real_t conicxx_max_step_fraction{(real_t)0.99};     // conicxx.max_step_fraction -- fraction-to-boundary safety factor
    bool conicxx_equilibrate{true};                     // conicxx.equilibrate -- enable Ruiz scaling
    int conicxx_equilibrate_max_iter{10};               // conicxx.equilibrate_max_iter
    real_t conicxx_equilibrate_min_scale{(real_t)1e-4};  // conicxx.equilibrate_min_scale
    real_t conicxx_equilibrate_max_scale{(real_t)1e4};   // conicxx.equilibrate_max_scale
    bool conicxx_validate_inputs{true};                 // conicxx.validate_inputs

    // Presence flags for precedence handling (e.g., pj.alpha overrides alpha if both set)
    bool has_pj_alpha{false};
    bool has_pj_max_iterations{false};
    bool has_pj_tol_abs{false};
    bool has_pj_tol_rel{false};
    bool has_pj_relaxation{false};

    // Integrator selection
    IntegratorType integrator{IntegratorType::Moreau};
    real_t moreau_theta{(real_t)1.0};
    bool moreau_implicit_gyroscopy{false};
    bool moreau_lambda_theta{false};

    std::string scene_name{"none-specified"};
};

class ConfigReader {
   public:
    // Load config from file. If loading/parsing fails, returns a Config with header defaults.
    static Config fromFile(const std::string& path);
};

}  // namespace cardillo::config
