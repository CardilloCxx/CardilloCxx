// Natural map of the contact law: the common measure of the final check (solver/natural_map.hpp),
// zero exactly at solutions of Coulomb's and of the convexified law, its values for known
// violations and invariances; the stopping test of the projection solvers with their own step
// (misc::naturalMapError2); projected Jacobi, projected Gauss-Seidel and the condensed solver
// stopping on it, checked by the final check.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <cardillo/config/config.hpp>
#include <cardillo/misc/contact_projection.hpp>
#include <cardillo/physics/api/physics_engine.hpp>
#include <cardillo/physics/assembly/dynamics_assembler.hpp>
#include <cardillo/physics/solver/natural_map.hpp>

using namespace cardillo;
using solver::ContactNaturalMap;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

bool near(real_t a, real_t b, real_t tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

CsrMatrix sparse(const MatrixXXr& dense) {
    return dense.sparseView();
}

VectorXr vec(std::initializer_list<real_t> values) {
    VectorXr v(values.size());
    int i = 0;
    for (real_t x : values) v[i++] = x;
    return v;
}

// One frictional contact with W = I and unit inverse mass (rho = 1), mu = 0.5.
constexpr real_t kMu = 0.5;
ContactNaturalMap unitContact() {
    return ContactNaturalMap(sparse(Matrix33r::Identity()), VectorXr::Ones(3), VectorXr::Constant(3, kMu), 0);
}

void testSolutionsAreZeros() {
    const ContactNaturalMap map = unitContact();
    check(map.residual(vec({1, 0.1, -0.2}), vec({0, 0, 0})) == 0, "stick: percussion inside the cone, xi = 0");
    check(near(map.residual(vec({1, -kMu, 0}), vec({0, 2, 0})), 0), "slide: friction opposite to the slip, on the cone boundary");
    check(map.residual(vec({0, 0, 0}), vec({1, 0.3, 0.2})) == 0, "open contact");

    // Convexified law with the shift mu ||gamma_T^-|| of the pre-step slip 1: sliding at slip 2 lifts
    // the contact off with xi_N = mu (2 - 1). Coulomb's law is violated by the same point.
    const VectorXr shift = vec({kMu * 1, 0, 0});
    const VectorXr lambda = vec({1, -kMu, 0}), xi = vec({kMu, 2, 0});
    check(near(map.residual(lambda, xi, shift), 0), "convexified law: zero at its own solution");
    check(map.residual(lambda, xi) > 0.1 * kMu, "convexified solution violates Coulomb's law");

    // Frictionless rows: closed and open.
    const ContactNaturalMap normal(sparse(MatrixXXr::Identity(2, 2)), VectorXr::Ones(2), VectorXr::Zero(2), 2);
    check(normal.residual(vec({1, 0}), vec({0, 1})) == 0, "frictionless: closed and open contact");
}

void testViolations() {
    const ContactNaturalMap map = unitContact();
    // Stick with a normal velocity error: the residual is the velocity error itself.
    check(near(map.residual(vec({1, 0.1, 0}), vec({1e-3, 0, 0})), 1e-3), "stick: residual = normal velocity error");
    check(near(map.residual(vec({1, 0.1, 0}), vec({0, 0, 2e-3})), 2e-3), "stick: tangential velocity error");

    // Spurious percussion on a separating frictionless contact: the velocity it would cause,
    // G lambda, independent of the masses (G = 2: lambda = 0.2 gives 0.4 m/s; G = 20: lambda = 0.02).
    for (real_t G : {2.0, 20.0}) {
        const ContactNaturalMap normal(sparse(MatrixXXr::Identity(1, 1)), VectorXr::Constant(1, G), VectorXr::Zero(1), 1);
        check(near(normal.residual(vec({0.4 / G}), vec({1})), 0.4), "separating contact: residual G lambda, independent of the mass");
    }

    check(std::isnan(map.residual(vec({1, 0, std::numeric_limits<real_t>::quiet_NaN()}), vec({0, 0, 0}))), "NaN input gives NaN");
}

// Rotating the tangent frame of a contact (rows of W, lambda_T, xi_T) leaves the residual unchanged,
// also for an anisotropic Delassus block (rho = 1 / lambda_max(G) is invariant).
void testTangentFrameInvariance() {
    Matrix33r W;
    W << 1, 0.2, -0.3, 0.1, 1.5, 0.4, -0.2, 0.3, 0.8;
    const VectorXr Minv = vec({1.0, 2.0, 0.5});
    const VectorXr lambda = vec({0.7, 0.4, -0.5}), xi = vec({0.05, 0.3, -0.1});

    constexpr real_t angle = 0.73;
    Matrix33r R = Matrix33r::Identity();
    R.block<2, 2>(1, 1) << std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle);

    const ContactNaturalMap map(sparse(W), Minv, VectorXr::Constant(3, kMu), 0);
    const ContactNaturalMap rotated(sparse(R * W), Minv, VectorXr::Constant(3, kMu), 0);
    const real_t a = map.residual(lambda, xi);
    const real_t b = rotated.residual(R * lambda, R * xi);
    check(a > 0 && near(a, b), "tangent-frame rotation leaves the residual unchanged");
}

// Stopping test of the projection solvers: the increment of one projection step divided by the
// step, r = (lambda - lambda^+) / rho, at the input (lambda, xi) of the step.
Vector3r stepResidual(const Vector3r& lambda, const Vector3r& xi, real_t rho) {
    Vector3r next;
    misc::proxCoulomb(rho * xi - lambda, kMu, next);
    return (lambda - next) / rho;
}

void testProjectionStep() {
    // Solutions: open contact, sliding with friction opposite to the slip on the cone boundary.
    check(stepResidual(Vector3r::Zero(), Vector3r(1, 0.3, 0.2), 0.1).norm() == 0, "step: open contact");
    check(near(stepResidual(Vector3r(1, -kMu, 0), Vector3r(0, 2, 0), 0.1).norm(), 0), "step: sliding contact");

    // Stick: the residual is the velocity error itself, for any step (PJ's alpha).
    const Vector3r lambda(1, 0.1, 0), xi(1e-3, 0, 2e-3);
    for (real_t rho : {0.1, 0.01}) check((stepResidual(lambda, xi, rho) - xi).norm() < 1e-12, "step: residual = velocity error in stick, rho " + std::to_string(rho));

    // Tolerances: a normal error of a fast sliding contact is not hidden behind its slip;
    // e = 1e-3 / (1e-5 + 1e-2 * 1e-3) = 50.
    const real_t e2 = misc::naturalMapError2(1e-3, 0, 0, 1e-3, 1, 0, 1e-5, 1e-2);
    check(near(std::sqrt(e2), 50, 1e-9), "error: normal residual scaled by the normal velocity, not the slip");
    check(near(misc::naturalMapError2(2e-3, -1.0, 1e-3, 1e-3), 1), "error: frictionless row");

    // Non-finite iterates count as infinitely far from convergence.
    const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
    check(std::isinf(misc::naturalMapError2(0, nan, 0, 0, 0, 0, 1, 0)), "error: NaN gives infinity");
    check(std::isinf(misc::naturalMapError2(nan, 0, 1, 0)), "error: NaN gives infinity (frictionless)");
}

config::Config ballConfig() {
    config::Config cfg;
    cfg.sim_dt = 1e-2;
    cfg.sim_T = 5 * cfg.sim_dt;
    cfg.sim_gravity = Vector3r(0, 0, -9.81);
    cfg.output_interval_steps = 0;
    cfg.friction_default_mu = 0.3;
    cfg.pj_max_iterations = 100000;
    return cfg;
}

// Ball sliding on a plane: PJ stops on the natural map with its own step at the input of its last
// sweep; the final check evaluates the common natural map at the result. Both are velocity errors
// and agree near convergence, also for a small step alpha, at which the velocity increment of a
// sweep is much smaller than the remaining error (velocity_rms at the same tolerance leaves up to
// 8e-8 m/s at alpha = 0.02).
void testProjectedJacobiNaturalMap(bool blockPreconditioner, real_t alpha, bool nesterov) {
    const std::string tag = std::string(blockPreconditioner ? "block" : "diagonal") + " preconditioner, alpha " + std::to_string(alpha) + (nesterov ? ", Nesterov" : ", standard loop");
    config::Config cfg = ballConfig();
    cfg.solver_criterion = config::TerminationCriterion::NaturalMap;
    cfg.solver_final_check = config::FinalCheck::NaturalMap;
    cfg.pj_rdiag_true_delassus = blockPreconditioner;
    cfg.pj_alpha = alpha;
    cfg.pj_nesterov = nesterov;
    cfg.pj_tol_abs = 1e-9;
    cfg.pj_tol_rel = 1e-9;

    physics::PhysicsEngine engine(cfg);
    engine.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    engine.addRigidBody(physics::SphereShape(0.5), physics::RigidState(Vector3r(0, 0, 0.5 - 1e-6), Vector3r(1, 0, 0)), physics::RigidProps(1));
    while (!engine.isFinished()) {
        engine.step();
        const auto& result = engine.lastStepResidual();
        check(result.has_value(), tag + ": final check evaluated");
        if (!result) return;
        check(!engine.dynamicsAssembler().contacts().empty() && engine.dynamicsAssembler().contacts().front().last_impulse[0] > 0, tag + ": active contact");
        // Tolerance of the sliding contact: 1e-9 (1 + |xi_T|) with |xi_T| <= 1 m/s.
        std::ostringstream value;
        value << result->contact;
        check(result->contact < 2e-9, tag + ": final natural map within the PJ tolerance (" + value.str() + " m/s)");
        check(result->contact == result->coulomb, tag + ": Coulomb's law is PJ's own law");
        check(result->bilateral == 0, tag + ": no bilateral rows");
    }
}

// The ball of testProjectedJacobiNaturalMap() with projected Gauss-Seidel and the condensed solver:
// they stop on the natural map of each block update (with the solver's step, or G_ii times the
// increment after a Newton local solve), at the state of that update.
void testBlockSolverNaturalMap(config::SolverType solver, bool nesterov, const std::string& sweepMode, const std::string& localSolve) {
    const std::string tag = std::string(solver == config::SolverType::ProjectedGaussSeidel ? "PGS" : "condensed " + sweepMode + " " + localSolve) + (nesterov ? ", Nesterov" : ", standard loop");
    config::Config cfg = ballConfig();
    cfg.solver = solver;
    cfg.solver_criterion = config::TerminationCriterion::NaturalMap;
    cfg.solver_final_check = config::FinalCheck::NaturalMap;
    cfg.pj_nesterov = nesterov;
    cfg.condensed_sweep_mode = sweepMode;
    cfg.condensed_local_solve = localSolve;
    cfg.pj_tol_abs = 1e-9;
    cfg.pj_tol_rel = 1e-9;

    physics::PhysicsEngine engine(cfg);
    engine.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    engine.addRigidBody(physics::SphereShape(0.5), physics::RigidState(Vector3r(0, 0, 0.5 - 1e-6), Vector3r(1, 0, 0)), physics::RigidProps(1));
    while (!engine.isFinished()) {
        engine.step();
        const auto& result = engine.lastStepResidual();
        check(result.has_value() && !engine.dynamicsAssembler().contacts().empty(), tag + ": active contact, final check evaluated");
        if (!result) return;
        std::ostringstream value;
        value << result->contact;
        check(result->contact < 2e-9, tag + ": final natural map within the tolerance (" + value.str() + " m/s)");
    }
}

// Double pendulum on perfect distance constraints, no contacts: projected Gauss-Seidel and the
// condensed solver without true_schur iterate the constraint rows, so the natural-map stop must
// cover them (their residual W^T u + b, scaled with W^T u + b), or it would stop after one sweep.
void testBilateralRowsNaturalMap(config::SolverType solver) {
    const std::string tag = solver == config::SolverType::ProjectedGaussSeidel ? "PGS double pendulum" : "condensed double pendulum";
    config::Config cfg = ballConfig();
    cfg.solver = solver;
    cfg.condensed_true_schur = false;
    cfg.collision_disable_all = true;
    cfg.solver_criterion = config::TerminationCriterion::NaturalMap;
    cfg.solver_final_check = config::FinalCheck::NaturalMap;
    cfg.pj_tol_abs = 1e-9;
    cfg.pj_tol_rel = 1e-9;

    physics::PhysicsEngine engine(cfg);
    const auto anchor = engine.addStaticBody(physics::SphereShape(0.01), physics::RigidState{});
    const auto first = engine.addPointMass(1.0, Vector3r(1, 0, 0), Vector3r::Zero(), 0.01);
    const auto second = engine.addPointMass(2.0, Vector3r(2, 0, 0), Vector3r::Zero(), 0.01);
    engine.addLinearDistanceConstraint(anchor, first, Vector3r::Zero(), Vector3r::Zero(), std::numeric_limits<real_t>::infinity(), 0.0, 1.0);
    engine.addLinearDistanceConstraint(first, second, Vector3r::Zero(), Vector3r::Zero(), std::numeric_limits<real_t>::infinity(), 0.0, 1.0);
    while (!engine.isFinished()) {
        engine.step();
        const auto& result = engine.lastStepResidual();
        check(result.has_value(), tag + ": final check evaluated");
        if (!result) return;
        std::ostringstream value;
        value << result->bilateral;
        check(result->bilateral < 1e-8, tag + ": constraint rows converged (" + value.str() + " m/s)");
    }
}

// Final check disabled: nothing is stored. Enabled on a pendulum (perfect distance constraint, no
// contacts): the bilateral rows are solved exactly by PJ.
void testFinalCheckOptions() {
    {
        config::Config cfg = ballConfig();
        physics::PhysicsEngine engine(cfg);
        engine.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
        engine.addRigidBody(physics::SphereShape(0.5), physics::RigidState(Vector3r(0, 0, 0.5 - 1e-6)), physics::RigidProps(1));
        engine.step();
        check(!engine.lastStepResidual().has_value(), "final check disabled by default");
    }

    // Pendulum on a perfect distance constraint: the bilateral row residual of PJ is round-off.
    config::Config cfg = ballConfig();
    cfg.solver_final_check = config::FinalCheck::NaturalMap;
    cfg.collision_disable_all = true;
    physics::PhysicsEngine engine(cfg);
    const auto anchor = engine.addStaticBody(physics::SphereShape(0.01), physics::RigidState{});
    const auto mass = engine.addPointMass(1.0, Vector3r(1, 0, 0), Vector3r::Zero(), 0.01);
    engine.addLinearDistanceConstraint(anchor, mass, Vector3r::Zero(), Vector3r::Zero(), std::numeric_limits<real_t>::infinity(), 0.0, 1.0);
    engine.step();
    const auto& result = engine.lastStepResidual();
    check(result.has_value() && result->contact == 0 && result->bilateral < 1e-12, "pendulum: bilateral rows exact, no contacts");
}

void testConfigParsing() {
    const auto path = std::filesystem::temp_directory_path() / "cardillo_natural_map_test.config";
    {
        std::ofstream configFile(path);
        configFile << "solver.criterion = NaTuRaL_MaP\n";
        configFile << "solver.final_check = Natural_Map\n";
        check((bool)configFile, "write temporary config");
    }

    const config::Config cfg = config::ConfigReader::fromFile(path.string());
    std::error_code error;
    std::filesystem::remove(path, error);
    check(cfg.solver_criterion == config::TerminationCriterion::NaturalMap, "solver.criterion, case-insensitive value");
    check(cfg.solver_final_check == config::FinalCheck::NaturalMap, "solver.final_check, case-insensitive value");
    check(config::Config{}.solver_criterion == config::TerminationCriterion::VelocityRms, "velocity_rms remains the default criterion");
    check(config::Config{}.solver_final_check == config::FinalCheck::None, "final check disabled by default");
}

}  // namespace

int main() {
    testSolutionsAreZeros();
    testViolations();
    testTangentFrameInvariance();
    testProjectionStep();
    testProjectedJacobiNaturalMap(false, 0.3, true);
    testProjectedJacobiNaturalMap(true, 0.3, true);
    testProjectedJacobiNaturalMap(false, 0.02, true);
    testProjectedJacobiNaturalMap(false, 0.02, false);
    for (bool nesterov : {true, false}) {
        testBlockSolverNaturalMap(config::SolverType::ProjectedGaussSeidel, nesterov, "gauss_seidel", "projection");
        testBlockSolverNaturalMap(config::SolverType::Condensed, nesterov, "gauss_seidel", "projection");
    }
    testBlockSolverNaturalMap(config::SolverType::Condensed, true, "colored", "projection");
    testBlockSolverNaturalMap(config::SolverType::Condensed, true, "gauss_seidel", "newton");
    testBilateralRowsNaturalMap(config::SolverType::ProjectedGaussSeidel);
    testBilateralRowsNaturalMap(config::SolverType::Condensed);
    testFinalCheckOptions();
    testConfigParsing();
    if (failures == 0) std::cout << "test_natural_map: all checks passed\n";
    return failures == 0 ? 0 : 1;
}
