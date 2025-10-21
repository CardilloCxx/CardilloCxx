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

    // Presence flags for precedence handling (e.g., pj.alpha overrides alpha if both set)
    bool has_pj_alpha{false};
    bool has_pj_max_iterations{false};
    bool has_pj_tol_abs{false};
    bool has_pj_relaxation{false};
    bool has_pj_compliance{false};

    // Simulation settings
    real_t sim_T{(real_t)5.0};   // sim.T - total simulation time
    real_t sim_dt{(real_t)1e-3}; // sim.dt - time step

    // Output settings
    int output_interval_steps{50}; // output.interval_steps - VTK write interval in steps
};

class ConfigReader {
public:
    // Load config from file. If loading/parsing fails, returns a Config with header defaults.
    static Config fromFile(const std::string& path);
};

}
