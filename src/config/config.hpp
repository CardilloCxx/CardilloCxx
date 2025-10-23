#pragma once

#include <string>
#include "../misc/types.hpp"

namespace cardillo::config {

struct Config {

    // PROJECTED_JACOBI scoped settings
    int    pj_max_iterations{200000};   // pj.max_iterations
    real_t pj_tol_abs{(real_t)1e-4};    // pj.tol_abs
    real_t pj_relaxation{(real_t)0.9};  // pj.relaxation
    real_t pj_alpha{(real_t)0.3};       // pj.alpha
    real_t pj_compliance{(real_t)1e-6}; // pj.compliance
    bool   pj_nesterov{false};          // pj.nesterov (enable Nesterov acceleration)
    real_t pj_nesterov_beta_threshold{(real_t)0.995}; // pj.nesterov_beta_threshold
    int    pj_nesterov_restart_limit{4};              // pj.nesterov_restart_limit
    bool   pj_warmstart{true};                        // pj.warmstart (enable warmstart & cache)

    // Presence flags for precedence handling (e.g., pj.alpha overrides alpha if both set)
    bool has_pj_alpha{false};
    bool has_pj_max_iterations{false};
    bool has_pj_tol_abs{false};
    bool has_pj_relaxation{false};
    bool has_pj_compliance{false};

    // Simulation settings
    real_t sim_T{(real_t)5.0};   // sim.T - total simulation time
    real_t sim_dt{(real_t)1e-3}; // sim.dt - time step
    // Gravity vector (world frame) used by examples and PhysicsSystem
    Vector3r sim_gravity{(real_t)0, (real_t)0, (real_t)-9.81}; // sim.gravity = gx gy gz

    // Output settings
    int output_interval_steps{50}; // output.interval_steps - VTK write interval in steps
    // Contact VTK detail: when true, include body-space normals/tangents/points in contacts VTK
    bool output_contacts_body_vectors{false}; // output.contacts_body_vectors

    // Collision settings
    // collision.broadphase: one of [dynamic_aabb, dynamic_aabb_array, naive, sap, ssap]
    std::string collision_broadphase{"dynamic_aabb"};
    // Maximum distance (in meters) to consider a current contact matched to a previous one
    // Used for warmstarting across steps. Set via `collision.match_max_dist`.
    real_t collision_match_max_dist{(real_t)0.02};
    // Friction settings
    bool friction_enable{false};         // friction.enable - extend W with tangents
    std::string friction_combine{"min"}; // friction.combine - one of [min, arith, geom]
    real_t friction_default_mu{(real_t)0.5}; // friction.default_mu - default C_Friction when not set on entities

    // Debug/diagnostics
    bool debug_rb{false};   // debug.rb  - enable rigid-body contact/W diagnostics in Moreau
    bool debug_pj{false};   // debug.pj  - enable ProjectedJacobi iteration logging
};

class ConfigReader {
public:
    // Load config from file. If loading/parsing fails, returns a Config with header defaults.
    static Config fromFile(const std::string& path);
};

}
