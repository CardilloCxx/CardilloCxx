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
    bool   pj_nesterov{false};          // pj.nesterov (enable Nesterov acceleration)
    real_t pj_nesterov_beta_threshold{(real_t)0.995}; // pj.nesterov_beta_threshold
    int    pj_nesterov_restart_limit{4};              // pj.nesterov_restart_limit
    bool   pj_warmstart{true};                        // pj.warmstart (enable warmstart & cache)

    // Presence flags for precedence handling (e.g., pj.alpha overrides alpha if both set)
    bool has_pj_alpha{false};
    bool has_pj_max_iterations{false};
    bool has_pj_tol_abs{false};
    bool has_pj_relaxation{false};
    

    // Simulation settings
    real_t sim_T{(real_t)5.0};   // sim.T - total simulation time
    real_t sim_dt{(real_t)1e-3}; // sim.dt - time step
    // Gravity vector (world frame) used by examples and PhysicsSystem
    Vector3r sim_gravity{(real_t)0, (real_t)0, (real_t)-9.81}; // sim.gravity = gx gy gz

    // Output settings
    int output_interval_steps{50}; // output.interval_steps - VTK write interval in steps
    int output_heightfield_stride{8}; // output.heightfield_stride - decimation factor for HeightField VTK
    bool output_contacts_body_vectors{false}; // output.contacts_body_vectors
    std::string output_folder{"./vtk_out"}; // output.folder
    std::string output_filename_prefix{"scene"}; // output.filename_prefix
    bool output_write_contacts{false}; // output.write_contacts

    // Collision settings
    // collision.broadphase: one of [dynamic_aabb, dynamic_aabb_array, naive, sap, ssap]
    std::string collision_broadphase{"dynamic_aabb"};
    // Disable all collisions globally (no collidable components, no contact computation)
    bool collision_disable_all{false};
    // Security margin for collision request (>=0). Small positive helps robustness.
    real_t collision_security_margin{(real_t)0.0}; // collision.security_margin
    // Maximum raw contacts requested per pair from narrowphase
    int collision_max_raw_contacts{1024};      // collision.max_raw_contacts
    // Maximum number of contact patches requested per pair
    std::size_t collision_max_patches{4};      // collision.max_patches
    // Prefer using patch vertices over raw contacts when available
    bool collision_use_patch_vertices{true};   // collision.use_patch_vertices
    // Maximum distance (in meters) to consider a current contact matched to a previous one
    // Used for warmstarting across steps. Set via `collision.match_max_dist`.
    real_t collision_match_max_dist{(real_t)0.02};
    // Minimum separation between contacts within the same entity-entity pair (meters).
    // Contacts closer than this will be deduplicated (keeping deeper penetrations).
    // Set via `collision.min_pair_contact_distance`. Set to 0 to disable.
    real_t collision_min_pair_contact_distance{(real_t)0.0};
    // Friction settings
    bool friction_enable{false};         // friction.enable - extend W with tangents
    std::string friction_combine{"min"}; // friction.combine - one of [min, arith, geom]
    real_t friction_default_mu{(real_t)0.5}; // friction.default_mu - default C_Friction when not set on entities

    // Debug/diagnostics
    bool debug_rb{false};   // debug.rb  - enable rigid-body contact/W diagnostics in Moreau
    bool debug_pj{false};   // debug.pj  - enable ProjectedJacobi iteration logging
    bool debug_mesh{false}; // debug.mesh - print mesh normalization info (volume, COM, inertia)

    std::string scene_name{"none-specified"};
};

class ConfigReader {
public:
    // Load config from file. If loading/parsing fails, returns a Config with header defaults.
    static Config fromFile(const std::string& path);
};

}
