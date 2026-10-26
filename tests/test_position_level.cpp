// Validation of the position-level force law against the theory in
// the paper (Section "A position-level variant without initial multipliers"):
//  1. rigid pendulum (distance constraint, C = 0): bounded constraint violation without any
//     Baumgarte parameter and without initial multipliers;
//  2. linear spring, theta = 1/2: the scheme is the trapezoidal rule, energy is conserved;
// All solvers are tested; they share the discretization of DynamicsAssembler::springBias().

#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <string>
#include <utility>

#include "physics/api/physics_engine.hpp"
#include "physics/assembly/dynamics_assembler.hpp"
#include "rigid_body/rigid_body.hpp"

using namespace cardillo;
using namespace cardillo::physics;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) {
        ++g_failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

struct Run {
    real_t maxViolation{0};  // max |g| over the run (distance constraint / spring elongation)
    real_t energyDrift{0};   // relative change of the total energy
    Vector3r finalPosition{Vector3r::Zero()};
    real_t percussion{0};    // bilateral percussion of the last step
};

// Point mass m = 1 attached to a static anchor at the origin by a distance law with rest length 1.
// stiffness = inf: rigid pendulum under gravity, released horizontally. Finite stiffness: spring
// along x without gravity, released with elongation 0.1.
Run simulateImpl(config::SolverType solver, real_t theta, real_t stiffness);

// Returns NaN results instead of aborting the whole test if a solver throws.
Run simulate(config::SolverType solver, real_t theta, real_t stiffness) {
    try {
        return simulateImpl(solver, theta, stiffness);
    } catch (const std::exception& ex) {
        if (std::string(ex.what()).find("Failed to load QOCO runtime backend") != std::string::npos) throw;  // environment, not the scheme
        std::cerr << "exception (theta=" << theta << ", stiffness=" << stiffness << "): " << ex.what() << std::endl;
        Run r;
        r.maxViolation = r.energyDrift = r.percussion = std::numeric_limits<real_t>::quiet_NaN();
        r.finalPosition.setConstant(std::numeric_limits<real_t>::quiet_NaN());
        return r;
    }
}

Run simulateImpl(config::SolverType solver, real_t theta, real_t stiffness) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.moreau_theta = theta;
    cfg.sim_dt = 1e-2;
    cfg.sim_T = 10.0;
    cfg.output_interval_steps = 0;
    cfg.collision_disable_all = true;
    cfg.pj_tol_abs = 1e-12;
    cfg.pj_tol_rel = 1e-12;
    cfg.pj_max_iterations = 100000;
    cfg.condensed_true_schur = true;  // pure DAE: exact bilateral solve
    const bool spring = std::isfinite(stiffness);
    cfg.sim_gravity = spring ? Vector3r::Zero() : Vector3r(0, 0, -9.81);

    PhysicsEngine engine(cfg);
    const auto anchor = engine.addStaticBody(SphereShape(0.01), RigidState{});
    const auto mass = engine.addPointMass(1.0, Vector3r(spring ? 1.1 : 1.0, 0, 0), Vector3r::Zero(), 0.01);
    const size_t idx = engine.addLinearDistanceConstraint(anchor, mass, Vector3r::Zero(), Vector3r::Zero(), stiffness, 0.0, 1.0);
    auto& reg = engine.ecs();
    auto& constraint = *engine.world().constraintPatterns()[idx];

    auto energy = [&]() {
        const Vector3r x = reg.get<C_Position3>(mass).value;
        const Vector3r v = reg.get<C_LinearVelocity3>(mass).value;
        const real_t e = x.norm() - 1.0;
        return 0.5 * v.squaredNorm() - cfg.sim_gravity.dot(x) + (spring ? 0.5 * stiffness * e * e : 0.0);
    };
    const real_t E0 = energy();
    const real_t Escale = spring ? E0 : 9.81;  // pendulum: relative to m g L

    Run r;
    while (!engine.isFinished()) {
        engine.step();
        RigidBody::updateState(reg, mass);  // constraints read the cached pose
        r.maxViolation = std::max(r.maxViolation, std::abs(constraint.getConstraint().positionError[0]));
    }
    r.energyDrift = (energy() - E0) / Escale;
    r.finalPosition = reg.get<C_Position3>(mass).value;
    r.percussion = engine.dynamicsAssembler().Lambda_g()[0];
    return r;
}

}  // namespace

int main() {
    using config::SolverType;
    const std::pair<SolverType, std::string> solvers[] = {{SolverType::ProjectedJacobi, "projected_jacobi"}, {SolverType::Condensed, "condensed"}, {SolverType::ProjectedGaussSeidel, "pgs"},
                                                          {SolverType::ConjugateGradient, "cg"},          {SolverType::Qoco, "qoco"},           {SolverType::Clarabel, "clarabel"},
                                                          {SolverType::Conicxx, "conicxx"}};
    // Bilateral percussions of projected Jacobi (first solver) as reference for all solvers, in
    // particular for the perfect constraint, whose percussion is a multiplier in the conic solvers.
    real_t refPercussion[3] = {0, 0, 0};
    for (const auto& [solver, sname] : solvers) {
        auto checkPercussion = [&](int idx, const Run& r, const std::string& what) {
            if (solver == SolverType::ProjectedJacobi) refPercussion[idx] = r.percussion;
            check(std::abs(r.percussion - refPercussion[idx]) < 1e-6 * std::abs(refPercussion[idx]) + 1e-9,
                  sname + ": " + what + " percussion " + std::to_string(r.percussion) + " as projected Jacobi " + std::to_string(refPercussion[idx]));
        };
        try {
            // 1) rigid pendulum
            for (int i = 0; i < 2; ++i) {
                const real_t theta = i == 0 ? 0.6 : 1.0;
                const Run r = simulate(solver, theta, INFINITY);
                std::cout << sname << " pendulum theta=" << theta << ": max|g| " << r.maxViolation << ", percussion " << r.percussion << std::endl;
                check(r.maxViolation < 2e-3, sname + ": pendulum violation bounded (theta=" + std::to_string(theta) + ")");
                checkPercussion(i, r, "pendulum (theta=" + std::to_string(theta) + ")");
            }

            // 2) linear spring, theta = 1/2: energy conservation (trapezoidal rule)
            {
                const Run r = simulate(solver, 0.5, 100.0);
                std::cout << sname << " spring theta=0.5: rel. energy change " << r.energyDrift << std::endl;
                check(std::abs(r.energyDrift) < 1e-6, sname + ": spring conserves energy at theta=1/2");
                checkPercussion(2, r, "spring");
            }
        } catch (const std::exception& ex) {
            std::cout << sname << ": SKIPPED (" << std::string(ex.what()).substr(0, 60) << ")" << std::endl;
        }
    }
    if (g_failures == 0) std::cout << "test_position_level: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
