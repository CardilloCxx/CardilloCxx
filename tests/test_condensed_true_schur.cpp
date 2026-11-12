// Exact elimination of the bilateral subsystem in the condensed solver (condensed.true_schur)
// against the direct solve of the projected Jacobi solver.
//
// The test uses a damped rod hanging from a static support: every rod element contributes a
// spring and a damper block between the same two nodes, i.e. pairs of bilateral blocks that share
// two bodies. Their coupling in Sbb must be the sum of the contributions of both bodies, each
// counted once (it was counted twice before, which made true_schur diverge, e.g. in the cable
// example).
//  1. Without contacts, both solvers solve the same linear one-step problem exactly, so their
//     trajectories must coincide up to round-off.
//  2. With contacts (the rod swings onto a static sphere) and a loose tolerance of the contact
//     iteration, the returned velocities must still satisfy the bilateral rows exactly for the
//     final percussions, as those of projected Jacobi do.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <cardillo/misc/spline.hpp>
#include <cardillo/physics/api/physics_engine.hpp>
#include <cardillo/physics/assembly/dynamics_assembler.hpp>

using namespace cardillo;

namespace {

struct Run {
    std::vector<Vector3r> positions;
    real_t maxBilateral{0};  // max |W^T u + b + c pi| over all spring and damper rows and steps
    int contactSteps{0};     // steps with contacts between the rod and the sphere
};

// Residual of the spring and damper rows after a step. The assembler still holds the rows, biases and
// percussions of that step; the bodies hold the new velocities u_{n+1}.
real_t bilateralResidual(physics::PhysicsEngine& engine) {
    const auto& dyn = engine.dynamicsAssembler();
    const auto& offsets = dyn.bodyVelOffsets();
    VectorXr u = VectorXr::Zero(offsets.back());
    for (auto [e, bi] : engine.ecs().view<C_BodyIndex, C_PhysicsObject>().each()) {
        if (bi.b < 0 || bi.b + 1 >= (int)offsets.size()) continue;
        const int nV = offsets[(size_t)bi.b + 1] - offsets[(size_t)bi.b];
        if (nV > 0) u.segment(offsets[(size_t)bi.b], nV) = engine.world().getVelocity(e).head(nV);
    }

    const real_t dt = engine.timeStep();
    const real_t theta = engine.config().moreau_theta;
    const real_t cSpring = physics::DynamicsAssembler::springComplianceScale(dt, theta);
    const real_t cDamper = physics::DynamicsAssembler::damperComplianceScale(dt, theta);
    real_t r = 0;
    if (dyn.Cdiag().size() > 0) r = std::max(r, (dyn.Wg().asSparse() * u + dyn.springBias(dt, theta) + cSpring * dyn.Cdiag().cwiseProduct(dyn.Lambda_g())).cwiseAbs().maxCoeff());
    if (dyn.Adiag().size() > 0) r = std::max(r, (dyn.Wgamma().asSparse() * u + dyn.damperBias(dt, theta) + cDamper * dyn.Adiag().cwiseProduct(dyn.Lambda_gamma())).cwiseAbs().maxCoeff());
    return r;
}

Run simulate(config::SolverType solver, bool trueSchur, bool obstacle) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.condensed_true_schur = trueSchur;
    cfg.moreau_theta = 0.6;
    cfg.sim_dt = 1e-3;
    cfg.sim_T = obstacle ? 0.3 : 0.2;
    cfg.output_interval_steps = 0;
    cfg.collision_disable_all = !obstacle;
    // Without contacts the one-step problem is linear and solved exactly; with contacts the loose
    // tolerance leaves the contact iteration far from converged.
    cfg.pj_tol_abs = obstacle ? 1e-3 : 1e-12;
    cfg.pj_tol_rel = obstacle ? 1e-3 : 1e-12;
    cfg.pj_max_iterations = 100000;

    physics::PhysicsEngine engine(cfg);
    const auto support = engine.addStaticBody(physics::SphereShape(0.01), physics::RigidState(Vector3r::Zero()));

    // Horizontal rod that swings down under gravity.
    misc::LinearSpline line(Vector3r::Zero(), Vector3r(0.5, 0, 0));
    const auto section = physics::BeamCrossSection::round(0.005);
    auto springs = physics::BeamSpringParams::fromMaterial(1e7, 0.3);
    springs.setDampingFromFactor(0.01);
    const physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());
    const auto ends = engine.createBeam(line, section, springs, stateDefaults, physics::RigidProps::withDensity(1000.0), 20);
    engine.addRigidConstraint(ends.first, support);
    engine.disableCollisionBetween(ends.first, support);

    // Sphere in the path of the rod, 25 degrees below the horizontal.
    const real_t angle = 25.0 * M_PI / 180.0;
    entt::entity sphere{entt::null};
    if (obstacle) sphere = engine.addStaticBody(physics::SphereShape(0.05), physics::RigidState(Vector3r(0.4 * std::cos(angle), 0, -0.4 * std::sin(angle))));

    Run r;
    while (!engine.isFinished()) {
        engine.step();
        r.maxBilateral = std::max(r.maxBilateral, bilateralResidual(engine));
        const auto& contacts = engine.dynamicsAssembler().contacts();
        if (std::any_of(contacts.begin(), contacts.end(), [&](const collision::Contact& c) { return c.a == sphere || c.b == sphere; })) ++r.contactSteps;
    }

    for (auto e : engine.ecs().view<C_BodyIndex>()) r.positions.push_back(engine.ecs().get<C_Position3>(e).value);
    return r;
}

real_t maxDifference(const std::vector<Vector3r>& a, const std::vector<Vector3r>& b) {
    if (a.size() != b.size()) return INFINITY;
    real_t d = 0;
    for (size_t i = 0; i < a.size(); ++i) d = std::max(d, (a[i] - b[i]).norm());
    return d;
}

}  // namespace

int main() {
    int failures = 0;
    auto check = [&](bool ok, const std::string& what) {
        if (!ok) {
            ++failures;
            std::cerr << "FAILED: " << what << std::endl;
        }
    };

    // 1) without contacts
    const Run pj = simulate(config::SolverType::ProjectedJacobi, false, false);
    const Run ts = simulate(config::SolverType::Condensed, true, false);
    real_t sag = 0;
    for (const auto& p : pj.positions) sag = std::max(sag, -p.z());
    const real_t diff = maxDifference(pj.positions, ts.positions);
    std::cout << "test_condensed_true_schur: rod tip sag " << sag << " m, max position difference condensed(true_schur) vs projected Jacobi " << diff << " m" << std::endl;
    check(sag > 0.05, "the rod does not swing down (sag " + std::to_string(sag) + ")");
    check(diff < 1e-8, "condensed with true_schur deviates from projected Jacobi by " + std::to_string(diff) + " m");

    // 2) with contacts: bilateral rows of the returned velocities, projected Jacobi as reference
    const Run pjc = simulate(config::SolverType::ProjectedJacobi, false, true);
    const Run tsc = simulate(config::SolverType::Condensed, true, true);
    std::cout << "test_condensed_true_schur: rod on the sphere in " << tsc.contactSteps << " steps, max bilateral residual condensed(true_schur) " << tsc.maxBilateral << " m/s, projected Jacobi "
              << pjc.maxBilateral << " m/s" << std::endl;
    check(tsc.contactSteps > 0, "the rod does not hit the sphere");
    check(pjc.maxBilateral < 1e-8, "projected Jacobi violates the bilateral rows by " + std::to_string(pjc.maxBilateral) + " m/s");
    check(tsc.maxBilateral < 1e-8, "condensed with true_schur violates the bilateral rows by " + std::to_string(tsc.maxBilateral) + " m/s");

    if (failures == 0) std::cout << "test_condensed_true_schur: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
