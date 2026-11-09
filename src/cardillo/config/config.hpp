#include <map>
#include <optional>
#pragma once

#include <string>
#include "cardillo/misc/types.hpp"

namespace cardillo::config {

enum class IntegratorType { Moreau };
// Discretization of the compliant force laws C*lambda + g(q) = 0 in the Moreau-theta scheme: the
// theta-weighted force lambda_theta := theta lambda_{n+1} + (1 - theta) lambda_{n} is the unknown, 
// and the position-level law is linearized at the intermediate configuration q_{n+theta} where 
// W is evaluated:
//   C*lambda_theta + g(t_{n+theta}, q_{n+theta}) + h*(theta^2 W^T u_{n+1} - (1-theta)^2 W^T u_n)
//     + (2 theta - 1)*h*chi = 0,
// with free bodies, kinematic drivers and prescribed rates all at t_{n+theta} = t_n + (1 - theta) h.
// No initial multipliers, no drift, no stabilization parameter; perfect constraints require
// theta > 1/2.
enum class SolverType { ProjectedJacobi, ProjectedGaussSeidel, Condensed, NaiveSchur, Qoco, Clarabel, Conicxx, Scs, Mosek };
enum class PjTerminationCriterion { VelocityRms, NaturalMap };
enum class PjFinalCheck { None, NaturalMap };

struct Config {
    // Simulation settings
    real_t sim_T{(real_t)5.0};    // sim.T - total simulation time
    real_t sim_dt{(real_t)1e-3};  // sim.dt - time step
    // Gravity vector (world frame) used by examples and World
    Vector3r sim_gravity{(real_t)0, (real_t)0, (real_t)-9.81};  // sim.gravity = gx gy gz

    // Output settings
    int output_interval_steps{50};                // output.interval_steps - VTK write interval in steps
    std::string output_folder{"./vtk_out"};       // output.folder
    std::string output_filename_prefix{""};  // output.filename_prefix (empty: outputPrefix() falls back to "scene")
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
    // Number of sampled support points for curved contact faces such as cylinder caps
    std::size_t collision_patch_num_samples{12};  // collision.patch_num_samples
    // Absolute distance from the support plane for including a point in a contact patch. Per pair,
    // it is raised to the penetration depth, so that every penetrating vertex is part of the patch.
    real_t collision_patch_tolerance{(real_t)1e-3};  // collision.patch_tolerance
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
    bool friction_enable{true};               // friction.enable - extend W with tangents
    std::string friction_combine{"arith"};    // friction.combine - one of [min, arith, geom]
    real_t friction_default_mu{(real_t)0.3};  // friction.default_mu - default C_Friction when not set on entities

    // Restitution settings
    real_t restitution_default_normal{(real_t)0.0};      // restitution.default_normal
    real_t restitution_default_tangential{(real_t)0.0};  // restitution.default_tangential

    // Debug/diagnostics
    bool debug_rb{false};    // debug.rb  - enable rigid-body contact/W diagnostics in Moreau
    bool debug_pj{false};    // debug.pj  - enable ProjectedJacobi iteration logging

    // Solver selection
    SolverType solver{SolverType::ProjectedJacobi};
    int pj_max_iterations{200000};                     // pj.max_iterations
    PjTerminationCriterion pj_criterion{PjTerminationCriterion::VelocityRms};  // pj.criterion
    PjFinalCheck pj_final_check{PjFinalCheck::None};    // pj.final_check
    real_t pj_tol_abs{(real_t)1e-5};                   // pj.tol_abs
    real_t pj_tol_rel{(real_t)1e-5};                   // pj.tol_rel
    // Step-length factors. PGS/condensed: dlambda = relaxation * GiiInv * r on every row, with an
    // additional factor alpha on contact rows. PJ ignores relaxation and uses the contact step
    // R = alpha * D^-1 (D = W M^-1 W^T); alpha bounds the step by the Lipschitz constant of the
    // contact Delassus operator. Both are tuned per scene.
    real_t pj_relaxation{(real_t)1.0};                 // pj.relaxation
    real_t pj_alpha{(real_t)0.3};                      // pj.alpha
    bool pj_nesterov{true};                            // pj.nesterov (enable Nesterov acceleration)
    real_t pj_nesterov_beta_threshold{(real_t)1.0};    // pj.nesterov_beta_threshold
    int pj_nesterov_restart_limit{10000};              // pj.nesterov_restart_limit
    // Chebyshev semi-iterative acceleration with a power-iteration estimate of the spectral radius
    // of the linearized sweep (PJ and condensed). Precedence: pj_nesterov > pj_chebyshev > pj_anderson.
    bool pj_chebyshev{false};  // pj.chebyshev
    // Anderson acceleration (type II) over the last pj_anderson_m residual differences (PJ and naive Schur).
    bool pj_anderson{false};  // pj.anderson
    int pj_anderson_m{2};     // pj.anderson_m (history window size, >= 1)
    bool pj_warmstart{true};  // pj.warmstart (enable warmstart & cache)
    // pj.rdiag_true_delassus: false (default) uses per-row step lengths alpha/D_ii (with a common
    // value alpha/mean(D_t1t1, D_t2t2) for both tangential rows of a frictional contact); true uses
    // the "split" step lengths (1/D_NN, 1/lambda_max of the tangential block) per frictional contact.
    bool pj_rdiag_true_delassus{false};                // pj.rdiag_true_delassus

    // Condensed (matrix-free block-relaxation) solver settings. Iteration/tolerance/relaxation/
    // warmstart/debug knobs are shared with PJ/PGS -- see pj_max_iterations, pj_tol_abs,
    // pj_tol_rel, pj_relaxation, pj_alpha, pj_warmstart, debug_pj above.
    std::string condensed_sweep_mode{"gauss_seidel"};    // condensed.sweep_mode: gauss_seidel | colored
    std::string condensed_local_solve{"projection"};     // condensed.local_solve: projection | newton
    int condensed_num_threads{0};                        // condensed.num_threads (0 = OpenMP default)
    int condensed_newton_max_iters{8};                   // condensed.newton_max_iters
    real_t condensed_newton_tol{(real_t)1e-10};          // condensed.newton_tol
    // Step lengths of the local Newton solve: "split" (normal and tangential treated separately) or
    // "full" (one rho from the largest eigenvalue of the 3x3 block); Siconos'
    // compute_rho_split_spectral_norm / compute_rho_spectral_norm.
    std::string condensed_newton_rho_strategy{"split"}; // condensed.newton_rho_strategy: split | full
    // How the step lengths (rhoN, rhoT, rhoT) of the plain projection update (GiiInv, used by every
    // sweep mode's default/fallback path) are derived from the local Delassus block of
    // ContactFrictional blocks (misc/contact_rho.hpp): "diagonal" (default) rhoN = 1/Gii(0,0),
    // rhoT = 2/(Gii(1,1)+Gii(2,2)); "split" rhoT = 1/lambda_max of the tangential 2x2 block; "full" a
    // single 1/lambda_max of the whole block. In all cases both tangential rows share one step
    // length and normal and tangential rows are not coupled -- anything else changes the fixed point
    // of the projection, i.e. the converged impulses violate Coulomb's law (friction not opposite to
    // the slip velocity). Spring/damper rows always use the per-row diagonal 1/(Gii(i,i)+c_i).
    std::string condensed_projection_rho_strategy{"diagonal"}; // condensed.projection_rho_strategy: diagonal | split | full
    // condensed.true_schur: eliminate the bilateral (spring+damper) rows exactly in every sweep with
    // a block-sparse LDLT (CondensedAssembler::buildBilateralFactorization()) instead of sweeping them.
    bool condensed_true_schur{false};                    // condensed.true_schur


    // Interior-point solver settings shared across QOCO, Clarabel and ConicXX
    real_t ip_kkt_static_reg{(real_t)1e-8};   // solver.kkt_static_reg
    real_t ip_kkt_dynamic_reg{(real_t)1e-8};  // solver.kkt_dynamic_reg
    real_t ip_tol_abs{(real_t)1e-8};          // ip.tol_abs
    real_t ip_tol_rel{(real_t)1e-8};          // ip.tol_rel
    int ip_iter_ref_iters{2};                 // solver.iter_ref_iters

    // MOSEK-only (solver.type = mosek, requires -DCARDILLO_WITH_MOSEK=ON). Negative values and
    // "default" keep MOSEK's own defaults.
    int mosek_num_threads{0};               // mosek.num_threads (0 = let MOSEK decide)
    real_t mosek_tol_rel_gap{(real_t)-1};   // mosek.tol_rel_gap (MSK_DPAR_INTPNT_CO_TOL_REL_GAP)
    real_t mosek_tol_pfeas{(real_t)-1};     // mosek.tol_pfeas (MSK_DPAR_INTPNT_CO_TOL_PFEAS)
    real_t mosek_tol_dfeas{(real_t)-1};     // mosek.tol_dfeas (MSK_DPAR_INTPNT_CO_TOL_DFEAS)
    int mosek_max_iterations{-1};           // mosek.max_iterations (MSK_IPAR_INTPNT_MAX_ITERATIONS)
    real_t mosek_max_time{(real_t)-1};      // mosek.max_time [s] (MSK_DPAR_OPTIMIZER_MAX_TIME)
    std::string mosek_presolve{"default"};  // mosek.presolve [default, on, off, free]
    bool mosek_scaling{true};               // mosek.scaling -- diagonal pre-scaling (see mosek::ConicSolver)

    // ConicXX-only: unlike QOCO/Clarabel it keeps its solver across steps and warm-starts from the
    // previous solution while the active contact set is unchanged (see
    // docs/chapters/solvers/interior_point.rst).
    bool conicxx_warm_start{true};               // conicxx.warm_start
    real_t conicxx_warm_lambda{(real_t)0.99};    // conicxx.warm_lambda -- weight of the previous solution in the warm start
    bool conicxx_warm_start_after_setup{true};   // conicxx.warm_start_after_setup -- warm start from velocities and matched impulses when the contact set changed
    bool conicxx_compliant_rows{true};           // conicxx.compliant_rows -- springs/dampers as compliant zero-cone rows (off: percussions as primal variables)
    int conicxx_verbose{0};                      // conicxx.verbose -- 0 silent, 1 summary, 2 iteration table per solve

    // ConicXX-only: remaining conicxx::Settings fields not already covered by
    // the shared ip_kkt_*/pj_tol_*/pj_max_iterations fields above.
    real_t conicxx_dynamic_reg_eps{(real_t)1e-13};      // conicxx.dynamic_reg_eps -- pivots with D*sign below this are replaced by +-solver.kkt_dynamic_reg
    real_t conicxx_refine_tol{(real_t)1e-13};           // conicxx.refine_tol -- relative residual target of iterative refinement
    real_t conicxx_max_step_fraction{(real_t)0.99};     // conicxx.max_step_fraction -- fraction-to-boundary safety factor
    bool conicxx_equilibrate{true};                     // conicxx.equilibrate -- enable Ruiz scaling
    int conicxx_equilibrate_max_iter{10};               // conicxx.equilibrate_max_iter
    real_t conicxx_equilibrate_min_scale{(real_t)1e-4};  // conicxx.equilibrate_min_scale
    real_t conicxx_equilibrate_max_scale{(real_t)1e4};   // conicxx.equilibrate_max_scale
    bool conicxx_validate_inputs{true};                 // conicxx.validate_inputs

    // SCS-only (solver.type = scs, requires -DCARDILLO_WITH_SCS=ON). Negative numeric values keep
    // SCS's own defaults (scs_set_default_settings). SCS is a first-order (ADMM) solver, so its
    // tolerances are not comparable one-to-one with the interior-point ones; negative eps values fall
    // back to pj.tol_abs / pj.tol_rel like the other conic backends.
    real_t scs_eps_abs{(real_t)-1};              // scs.eps_abs
    real_t scs_eps_rel{(real_t)-1};              // scs.eps_rel
    real_t scs_eps_infeas{(real_t)-1};           // scs.eps_infeas
    int scs_max_iters{-1};                       // scs.max_iters
    real_t scs_time_limit_secs{(real_t)-1};      // scs.time_limit_secs
    real_t scs_alpha{(real_t)-1};                // scs.alpha (Douglas-Rachford relaxation)
    real_t scs_scale{(real_t)-1};                // scs.scale (initial dual scale)
    real_t scs_rho_x{(real_t)-1};                // scs.rho_x
    int scs_acceleration_lookback{-1};           // scs.acceleration_lookback (0 disables Anderson acceleration)
    int scs_acceleration_interval{-1};           // scs.acceleration_interval
    bool scs_adaptive_scale{true};               // scs.adaptive_scale
    bool scs_normalize{true};                    // scs.normalize (data equilibration)
    bool scs_warm_start{true};                   // scs.warm_start -- x from v_old/Lambda, y from the tracked contact impulses
    bool scs_carry_scale{true};                  // scs.carry_scale -- start from the previous step's final adapted scale
    std::string scs_stats_csv{""};               // scs.stats_csv -- per-step SCS statistics (empty = off)
    std::string scs_dump_dir{""};                // scs.dump_dir -- write each step's cone program + warm start (see tests/conic_replay.cpp)
    int scs_dump_every{1};                       // scs.dump_every -- dump every N-th step

    // Integrator selection
    IntegratorType integrator{IntegratorType::Moreau};
    real_t moreau_theta{(real_t)1.0};
    bool moreau_implicit_gyroscopy{false};

    std::string scene_name{"none-specified"};
    // Scene-specific parameters: every key "scene.<name> = value" except scene.name (e.g.
    // scene.waterfall_height = 0.05), read by the scenes via sceneParam().
    std::map<std::string, std::string> scene_params;

    // Value of the scene parameter `name` (key "scene.<name>"), or `fallback` if it is not set.
    // Throws if the value is not a number.
    real_t sceneParam(const std::string& name, real_t fallback) const;

    // Prefix of all output file names.
    std::string outputPrefix() const { return output_filename_prefix.empty() ? std::string("scene") : output_filename_prefix; }
};

class ConfigReader {
   public:
    // Loads a config file over the defaults of Config. Throws with file:line on a missing file or an
    // invalid value; unknown keys only warn.
    static Config fromFile(const std::string& path);
};

}  // namespace cardillo::config
