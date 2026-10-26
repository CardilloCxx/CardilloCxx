// Error of the convexified Coulomb law of the conic solvers (QOCO, Clarabel, ConicXX) for a solid
// sphere on a horizontal plane, against closed-form solutions of the discrete scheme (theory manual,
// section "Second-order cone formulation of the one-step problem"). In every step with an active
// contact, the conic solvers satisfy
//     xi_N = mu (||gamma_T^+|| - ||gamma_T^-||)                                             (*)
// instead of xi_N = 0 of the exact Coulomb law (projected Jacobi). With kappa = 1 + m rho^2 / Theta
// (7/2 for the lever arm rho = R) and a = kappa mu^2 / (1 + kappa mu^2):
//  1. Persistent sliding (ball placed on the plane with velocity v_0, no spin): the ball sinks with
//     the velocity -kappa mu^2 g h (1 - a^k), and the final penetration is mu h v_0. Both errors are
//     O(h); projected Jacobi keeps v_z = 0.
//  2. Impact of a spinning ball with restitution e: the rebound velocity is
//     ((e - kappa mu^2) V - kappa mu^2 g h) / (1 + kappa mu^2) for a sliding and e V - mu ||gamma_T^-||
//     for a sticking impact instead of e V. This error is O(1) in h.
// The contact point of Coal lies in the middle of the penetration, so the lever arm rho of the
// contact is slightly smaller than R; the closed forms use the lever arm of the engine's contact.

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "physics/api/physics_engine.hpp"
#include "physics/assembly/dynamics_assembler.hpp"

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

constexpr real_t kMass = 1.0;
constexpr real_t kRadius = 0.1;
constexpr real_t kGravity = 10.0;
constexpr real_t kMu = 0.3;
// Contacts whose normal percussion is below this threshold are separated; the interior-point methods
// return small positive percussions for them instead of zero.
constexpr real_t kActive = 1e-6;

struct Step {
    Vector3r r{Vector3r::Zero()}, v{Vector3r::Zero()}, omega{Vector3r::Zero()};  // state at the beginning of the step
    Vector3r vNew{Vector3r::Zero()}, omegaNew{Vector3r::Zero()};                 // velocities at the end of the step
    bool active{false};                                                          // contact with percussion in this step
    Vector3r point{Vector3r::Zero()};                                            // contact point of the step

    // Velocity of the material point of the ball at the contact point. The spin is about the world
    // y-axis, so body-fixed and spatial angular velocities coincide.
    Vector3r pointVelocity(bool post) const { return post ? Vector3r(vNew + omegaNew.cross(point - r)) : Vector3r(v + omega.cross(point - r)); }
    real_t gammaN(bool post) const { return pointVelocity(post).z(); }
    real_t slip(bool post) const { return pointVelocity(post).head<2>().norm(); }
    // 1 + m rho^2 / Theta with the lever arm rho of the contact point, Theta = 2/5 m R^2.
    real_t kappa() const { return 1 + (point - r).squaredNorm() / (0.4 * kRadius * kRadius); }
};

struct Run {
    std::vector<Step> steps;
    real_t finalHeight{0};
};

// Sphere above the plane z = 0 with center height z0, velocity v0 along x and spin omega0 about y.
Run simulate(config::SolverType solver, real_t dt, real_t T, real_t z0, real_t v0, real_t omega0, real_t e) {
    config::Config cfg;
    cfg.solver = solver;
    cfg.moreau_theta = 1.0;
    cfg.sim_dt = dt;
    cfg.sim_T = T;
    cfg.sim_gravity = Vector3r(0, 0, -kGravity);
    cfg.output_interval_steps = 0;
    cfg.friction_default_mu = kMu;
    cfg.restitution_default_normal = e;
    cfg.restitution_default_tangential = 0.0;
    // Projected Jacobi to round-off; the interior-point methods stall for much tighter tolerances.
    const bool conic = solver != config::SolverType::ProjectedJacobi;
    cfg.pj_tol_abs = conic ? 1e-10 : 1e-12;
    cfg.pj_tol_rel = conic ? 1e-10 : 1e-12;
    cfg.pj_max_iterations = 100000;

    PhysicsEngine engine(cfg);
    engine.addStaticBody(PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), RigidState{});
    const auto ball = engine.addRigidBody(SphereShape(kRadius), RigidState(Vector3r(0, 0, z0), Vector3r(v0, 0, 0), Vector3r(0, omega0, 0)), RigidProps(kMass));
    auto& reg = engine.ecs();

    Run run;
    while (!engine.isFinished()) {
        Step s;
        s.r = reg.get<C_Position3>(ball).value;
        s.v = reg.get<C_LinearVelocity3>(ball).value;
        s.omega = reg.get<C_AngularVelocity3>(ball).value;
        engine.step();
        s.vNew = reg.get<C_LinearVelocity3>(ball).value;
        s.omegaNew = reg.get<C_AngularVelocity3>(ball).value;
        const auto& contacts = engine.dynamicsAssembler().contacts();
        if (!contacts.empty() && contacts.front().last_impulse(0) > kActive) {
            s.active = true;
            s.point = contacts.front().point;
        }
        run.steps.push_back(s);
    }
    run.finalHeight = reg.get<C_Position3>(ball).value.z() - kRadius;
    return run;
}

const char* name(config::SolverType s) {
    switch (s) {
        case config::SolverType::ProjectedJacobi:
            return "pj";
        case config::SolverType::Qoco:
            return "qoco";
        case config::SolverType::Clarabel:
            return "clarabel";
        default:
            return "conicxx";
    }
}

// Largest residual of the exact law xi_N = 0 and of the convexified law (*) over all active steps.
struct Residuals {
    real_t exact{0}, convexified{0};
};
Residuals lawResiduals(const Run& run, real_t e) {
    Residuals res;
    for (const Step& s : run.steps) {
        if (!s.active) continue;
        const real_t xiN = s.gammaN(true) + e * s.gammaN(false);
        res.exact = std::max(res.exact, std::abs(xiN));
        res.convexified = std::max(res.convexified, std::abs(xiN - kMu * (s.slip(true) - s.slip(false))));
    }
    return res;
}

// 1. Ball placed on the plane (1e-8 deep, so that the contact is detected in the first step), sliding
// with v0 = 2 m/s without spin. The sliding phase ends at t = 2 v0 / (7 mu g) = 0.19 s.
void slidingBall(config::SolverType solver, real_t dt) {
    const real_t v0 = 2.0;
    const Run run = simulate(solver, dt, 0.4, kRadius - 1e-8, v0, 0.0, 0.0);
    const Residuals res = lawResiduals(run, 0.0);
    const std::string tag = std::string(name(solver)) + " sliding, h=" + std::to_string(dt);

    // Normal velocity after step k (1-based) in the middle of the sliding phase.
    const int k = static_cast<int>(std::lround(0.1 / dt));
    const Step& s = run.steps[k - 1];
    const real_t kmu2 = s.kappa() * kMu * kMu;
    const real_t vzConvex = -kmu2 * kGravity * dt * (1 - std::pow(kmu2 / (1 + kmu2), k));
    const real_t depth = -run.finalHeight;
    const real_t vxEnd = run.steps.back().vNew.x();
    std::cout << tag << ": max|xi_N| " << res.exact << ", max|xi_N - mu dslip| " << res.convexified << ", v_z(0.1 s) " << s.vNew.z() << " (convexified " << vzConvex << "), penetration " << depth
              << " (mu h v0 = " << kMu * dt * v0 << "), v_x(end) " << vxEnd << " (5/7 v0 = " << 5.0 / 7.0 * v0 << ")" << std::endl;

    if (solver != config::SolverType::ProjectedJacobi) {
        check(res.convexified < 1e-7, tag + ": convexified law (*) holds in every step");
        check(std::abs(s.vNew.z() - vzConvex) < 1e-3 * std::abs(vzConvex), tag + ": sinking velocity -kappa mu^2 g h (1 - a^k)");
        check(std::abs(depth - kMu * dt * v0) < 1e-2 * kMu * dt * v0, tag + ": final penetration mu h v0");
        check(std::abs(vxEnd - 5.0 / 7.0 * v0) < depth / kRadius * v0, tag + ": rolling velocity 5/7 v0 up to the lever arm R - depth/2");
    } else {
        check(res.exact < 1e-9, tag + ": exact Signorini law, v_z = 0");
        check(depth < 1e-7, tag + ": no penetration");
        check(std::abs(vxEnd - 5.0 / 7.0 * v0) < 1e-7, tag + ": rolling velocity 5/7 v0");
    }
}

// 2. Ball dropped from 0.2 m onto the plane with spin omega0 and restitution e: rebound velocity of
// the first impact. For e = 0, the ball stays in contact until it rolls (t < 1 s), and the convexified
// law (*) telescopes to the penetration mu h ||gamma_T^-|| after the impact (theta = 1).
void spinningImpact(config::SolverType solver, real_t dt, real_t omega0, real_t e) {
    const Run run = simulate(solver, dt, 1.0, kRadius + 0.2, 0.0, omega0, e);
    const std::string tag = std::string(name(solver)) + " impact, omega0=" + std::to_string(omega0) + " e=" + std::to_string(e) + " h=" + std::to_string(dt);

    const Step* impact = nullptr;
    for (const Step& s : run.steps)
        if (s.active) {
            impact = &s;
            break;
        }
    check(impact != nullptr, tag + ": impact detected");
    if (!impact) return;

    const real_t V = -impact->gammaN(false);
    const real_t slip = impact->slip(false);
    const real_t vPlus = impact->gammaN(true);
    const real_t kappa = impact->kappa();
    const real_t kmu2 = kappa * kMu * kMu;
    // Exact Coulomb impact: sticking if the slip can be removed within the friction cone.
    const bool stickExact = slip <= kappa * kMu * ((1 + e) * V + kGravity * dt);
    // Convexified impact: both branches meet at the transition, the sticking one has the larger rebound.
    const real_t vSlideConvex = ((e - kmu2) * V - kmu2 * kGravity * dt) / (1 + kmu2);
    const real_t vStickConvex = e * V - kMu * slip;
    const bool stickConvex = vStickConvex >= vSlideConvex;
    const real_t vConvex = stickConvex ? vStickConvex : vSlideConvex;
    std::cout << tag << ": V " << V << ", slip " << slip << ", kappa " << kappa << ", v+/V " << vPlus / V << " (exact " << e << ", " << (stickExact ? "stick" : "slide") << "; convexified "
              << vConvex / V << ", " << (stickConvex ? "stick" : "slide") << ")" << std::endl;

    // Penetration accumulated from the beginning of the impact step to the end of the run.
    const real_t depth = (impact->r.z() - kRadius) - run.finalHeight;
    if (e == 0) std::cout << tag << ": penetration after the impact " << depth << " (convexified mu h slip = " << kMu * dt * slip << ")" << std::endl;

    const Residuals res = lawResiduals(run, e);
    if (solver != config::SolverType::ProjectedJacobi) {
        check(res.convexified < 1e-6, tag + ": convexified law (*) holds in every step");
        check(std::abs(vPlus - vConvex) < 1e-6 * V, tag + ": rebound velocity of the convexified law");
        if (e == 0) check(std::abs(depth - kMu * dt * slip) < 2e-2 * kMu * dt * slip, tag + ": penetration mu h slip");
    } else {
        check(res.exact < 1e-9, tag + ": exact Signorini law");
        check(std::abs(vPlus - e * V) < 1e-9, tag + ": rebound velocity e V");
        if (e == 0) check(std::abs(depth) < 1e-9, tag + ": no penetration after the impact");
    }
}

}  // namespace

int main() {
    using config::SolverType;
    for (SolverType solver : {SolverType::ProjectedJacobi, SolverType::Qoco, SolverType::Clarabel, SolverType::Conicxx}) {
        try {
            for (real_t dt : {1e-3, 5e-4, 2.5e-4}) slidingBall(solver, dt);
            for (real_t dt : {1e-3, 1e-4}) {
                spinningImpact(solver, dt, 50.0, 0.5);  // sliding impact
                spinningImpact(solver, dt, 50.0, 0.0);
                spinningImpact(solver, dt, 10.0, 0.5);  // sticking impact
            }
        } catch (const std::exception& ex) {
            if (std::string(ex.what()).find("Failed to load QOCO runtime backend") != std::string::npos) {
                std::cout << name(solver) << ": skipped (" << ex.what() << ")" << std::endl;
                continue;
            }
            check(false, std::string(name(solver)) + ": exception " + ex.what());
        }
    }
    if (g_failures == 0) std::cout << "test_convexified_friction: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
