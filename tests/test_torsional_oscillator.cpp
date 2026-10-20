// Forced torsional oscillator on a revolute joint (examples/scenes/torsional_oscillator) against
// its closed-form solution. The bar rotates by more than 180 degrees, so the test checks that
//  - the hinge spring is linear in the rotation angle beyond 180 degrees (unwrapped twist angle of
//    the joint x-row),
//  - the position-level force law keeps the locked rows of the revolute joint satisfied for large
//    rotations, and
//  - the Moreau-theta scheme converges to the closed-form solution: first order with an error
//    proportional to (theta - 1/2) h, and joint violations of order h^2.
// theta = 1/2 is not tested: perfect constraints (here the locked joint rows, which carry the
// centripetal force) are not damped by the position-level law for theta = 1/2 and blow up.

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "../examples/scenes/torsional_oscillator/TorsionalOscillatorScene.hpp"
#include "physics/constraints/constraints.hpp"
#include "rigid_body/rigid_body.hpp"

using namespace cardillo;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) {
        ++g_failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

struct Run {
    real_t maxError{0};      // max |phi - phi_exact| over the run (rad)
    real_t maxAngle{0};      // max |phi| over the run (rad)
    real_t maxViolation{0};  // max violation of the locked joint rows
};

Run simulate(config::SolverType solver, real_t theta, real_t dt, real_t T) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.moreau_theta = theta;
    cfg.sim_dt = dt;
    cfg.sim_T = T;
    cfg.sim_gravity = Vector3r::Zero();
    cfg.output_interval_steps = 0;
    cfg.collision_disable_all = true;
    cfg.pj_tol_abs = 1e-12;
    cfg.pj_tol_rel = 1e-12;
    cfg.pj_max_iterations = 100000;
    cfg.condensed_true_schur = true;

    physics::PhysicsEngine engine(cfg);
    TorsionalOscillatorScene scene;
    scene.populate(engine);
    auto& reg = engine.ecs();
    const auto& joint = *engine.world().constraintPatterns().back();

    Run r;
    real_t t = 0;
    real_t phi = 0;  // unwrapped angle
    while (!engine.isFinished()) {
        scene.updateScene(engine, t);
        engine.step();
        t += dt;

        phi += std::remainder(scene.wrappedAngle(engine) - phi, (real_t)(2 * M_PI));
        r.maxError = std::max(r.maxError, std::abs(phi - scene.closedForm(t)));
        r.maxAngle = std::max(r.maxAngle, std::abs(phi));

        RigidBody::updateState(reg, scene.bar());  // constraints read the cached pose
        const VectorXr g = joint.getConstraint().positionError;
        for (int i : {0, 1, 2, 4, 5}) r.maxViolation = std::max(r.maxViolation, std::abs(g[i]));
    }
    return r;
}

}  // namespace

int main() {
    using config::SolverType;
    const real_t T = 5.0;
    // The error scales with the forcing amplitude A (static deflection M0 / k), the joint violation
    // (centripetal load) with A^2.
    const real_t A = TorsionalOscillatorScene().staticAngle;
    for (real_t theta : {0.55, 0.6}) {
        real_t prevError = 0;
        for (real_t dt : {2e-3, 1e-3, 5e-4}) {
            const Run pj = simulate(SolverType::ProjectedJacobi, theta, dt, T);
            const Run cond = simulate(SolverType::Condensed, theta, dt, T);
            const std::string tag = "theta=" + std::to_string(theta) + " dt=" + std::to_string(dt);
            std::cout << tag << ": max|phi - phi_exact| " << pj.maxError << " rad (condensed " << cond.maxError << "), max|phi| " << pj.maxAngle << " rad, max joint violation "
                      << pj.maxViolation << std::endl;

            check(pj.maxAngle > 1.5 * M_PI, tag + ": bar rotates by more than 270 degrees");
            check(pj.maxError < 300 * A * (theta - 0.5) * dt, tag + ": error of order (theta - 1/2) h");
            check(pj.maxViolation < 4 * A * A * dt * dt, tag + ": joint violation of order h^2");
            check(std::abs(pj.maxError - cond.maxError) < 1e-6 && std::abs(pj.maxViolation - cond.maxViolation) < 1e-8, tag + ": projected Jacobi and condensed agree");
            if (prevError > 0) {
                const real_t ratio = prevError / pj.maxError;
                check(ratio > 1.7 && ratio < 2.6, tag + ": first-order convergence (error ratio " + std::to_string(ratio) + ")");
            }
            prevError = pj.maxError;
        }
    }
    if (g_failures == 0) std::cout << "test_torsional_oscillator: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
