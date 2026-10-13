// Exact elimination of the bilateral subsystem in the condensed solver (condensed.true_schur)
// against the direct solve of the projected Jacobi solver. Without contacts, both solve the same
// linear one-step problem exactly, so their trajectories must coincide up to round-off.
//
// The test uses a damped rod hanging from a static support: every rod element contributes a
// spring and a damper block between the same two nodes, i.e. pairs of bilateral blocks that share
// two bodies. Their coupling in Sbb must be the sum of the contributions of both bodies, each
// counted once (it was counted twice before, which made true_schur diverge, e.g. in the cable
// example).

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "misc/spline.hpp"
#include "physics/api/physics_engine.hpp"

using namespace cardillo;

namespace {

std::vector<Vector3r> simulate(config::SolverType solver, bool trueSchur) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.condensed_true_schur = trueSchur;
    cfg.moreau_theta = 0.6;
    cfg.sim_dt = 1e-3;
    cfg.sim_T = 0.2;
    cfg.output_interval_steps = 0;
    cfg.collision_disable_all = true;
    cfg.pj_tol_abs = 1e-12;
    cfg.pj_tol_rel = 1e-12;
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

    while (!engine.isFinished()) engine.step();

    std::vector<Vector3r> positions;
    for (auto e : engine.ecs().view<C_BodyIndex>()) positions.push_back(engine.ecs().get<C_Position3>(e).value);
    return positions;
}

real_t maxDifference(const std::vector<Vector3r>& a, const std::vector<Vector3r>& b) {
    if (a.size() != b.size()) return INFINITY;
    real_t d = 0;
    for (size_t i = 0; i < a.size(); ++i) d = std::max(d, (a[i] - b[i]).norm());
    return d;
}

}  // namespace

int main() {
    const auto pj = simulate(config::SolverType::ProjectedJacobi, false);
    const auto ts = simulate(config::SolverType::Condensed, true);

    real_t sag = 0;
    for (const auto& p : pj) sag = std::max(sag, -p.z());
    const real_t diff = maxDifference(pj, ts);
    std::cout << "test_condensed_true_schur: rod tip sag " << sag << " m, max position difference condensed(true_schur) vs projected Jacobi " << diff << " m" << std::endl;

    int failures = 0;
    if (!(sag > 0.05)) {
        ++failures;
        std::cerr << "FAILED: the rod does not swing down (sag " << sag << ")" << std::endl;
    }
    if (!(diff < 1e-8)) {
        ++failures;
        std::cerr << "FAILED: condensed with true_schur deviates from projected Jacobi by " << diff << " m" << std::endl;
    }
    if (failures == 0) std::cout << "test_condensed_true_schur: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
