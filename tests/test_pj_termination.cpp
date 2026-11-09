#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <cardillo/config/config.hpp>
#include <cardillo/physics/assembly/dynamics_assembler.hpp>
#include <cardillo/physics/api/physics_engine.hpp>
#include <cardillo/physics/solver/termination.hpp>

using namespace cardillo;

namespace {

int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

void testNaturalMapResidual() {
    VectorXr lambdaOld(3), lambdaNew(3), stepSize(3);
    lambdaOld << 2, -1, 0.5;
    lambdaNew.setZero();
    stepSize << 0.5, 2, 0.25;

    const real_t error = solver::termination::naturalMapError(lambdaOld, lambdaNew, stepSize, VectorXr::Zero(3), VectorXr::Zero(3), 0, 1, 0);
    check(std::abs(error - 4) < 1e-12, "natural-map error applies each row's PJ step size");
}

void testBlockScaling() {
    VectorXr lambdaOld(4), lambdaNew(4), stepSize(4), rowVelocity(4), previousRowVelocity(4);
    lambdaOld << 0.2, 0.3, 0.4, 0;
    lambdaNew.setZero();
    stepSize.setOnes();
    rowVelocity << 0, 0, 1, 1;
    previousRowVelocity.setZero();

    const real_t error = solver::termination::naturalMapError(lambdaOld, lambdaNew, stepSize, rowVelocity, previousRowVelocity, 1, 0.1, 0.1);
    check(std::abs(error - 3) < 1e-12, "block error takes the maximum of normal and tangent errors");

    VectorXr contactOld(3), contactNew(3), contactStep(3), slidingVelocity(3);
    contactOld << 1, 0, 0;
    contactNew.setZero();
    contactStep.setOnes();
    slidingVelocity << 0, 1, 0;
    const real_t normalError = solver::termination::naturalMapError(contactOld, contactNew, contactStep, slidingVelocity, VectorXr::Zero(3), 0, 0.1, 0.1);
    check(std::abs(normalError - 10) < 1e-12, "sliding speed does not dilute normal error");
}

void testTangentFrameInvariance() {
    VectorXr lambdaOld(3), lambdaNew(3), stepSize(3), rowVelocity(3);
    lambdaOld << 0.1, 0.3, -0.4;
    lambdaNew.setZero();
    stepSize.setOnes();
    rowVelocity << 0.2, 1, 0.5;
    const real_t reference = solver::termination::naturalMapError(lambdaOld, lambdaNew, stepSize, rowVelocity, VectorXr::Zero(3), 0, 0.01, 0.1);

    constexpr real_t angle = 0.73;
    const real_t cosine = std::cos(angle);
    const real_t sine = std::sin(angle);
    const auto rotate = [cosine, sine](const VectorXr& value) {
        VectorXr rotated = value;
        rotated[1] = cosine * value[1] - sine * value[2];
        rotated[2] = sine * value[1] + cosine * value[2];
        return rotated;
    };
    const real_t rotatedError = solver::termination::naturalMapError(rotate(lambdaOld), lambdaNew, stepSize, rotate(rowVelocity), VectorXr::Zero(3), 0, 0.01, 0.1);
    check(std::abs(reference - rotatedError) < 1e-12, "tangent-frame rotation leaves block error unchanged");
}

void testNaturalMapContactSolve(bool blockPreconditioner) {
    config::Config cfg;
    cfg.pj_criterion = config::PjTerminationCriterion::NaturalMap;
    cfg.pj_final_check = config::PjFinalCheck::NaturalMap;
    cfg.debug_pj = true;
    cfg.pj_rdiag_true_delassus = blockPreconditioner;
    cfg.pj_tol_abs = 1e-10;
    cfg.pj_tol_rel = 1e-10;
    cfg.pj_max_iterations = 1000;
    cfg.sim_dt = 1e-2;
    cfg.sim_T = cfg.sim_dt;
    cfg.sim_gravity = Vector3r(0, 0, -9.81);
    cfg.output_interval_steps = 0;

    physics::PhysicsEngine engine(cfg);
    engine.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), 5, 5), physics::RigidState{});
    const auto ball = engine.addRigidBody(physics::SphereShape(0.5), physics::RigidState(Vector3r(0, 0, 0.49)), physics::RigidProps(1));
    std::ostringstream output;
    std::streambuf* previousOutput = std::cout.rdbuf(output.rdbuf());
    engine.step();
    std::cout.rdbuf(previousOutput);

    const auto& contacts = engine.dynamicsAssembler().contacts();
    check(!contacts.empty(), "natural-map PJ contact test generated a contact");
    if (!contacts.empty()) check(contacts.front().last_impulse[0] > 0, "natural-map PJ solve produced a positive normal impulse");
    check(engine.getLinearVelocity(ball).allFinite(), "natural-map PJ contact solve returned finite velocity");
    const std::string reportPrefix = "Final natural-map error: ";
    const std::size_t reportPos = output.str().find(reportPrefix);
    check(reportPos != std::string::npos, "enabled final check reports separately from the iteration error");
    if (reportPos != std::string::npos) {
        real_t finalError = 0;
        std::istringstream report(output.str().substr(reportPos + reportPrefix.size()));
        report >> finalError;
        check(report && std::isfinite(finalError), "final natural-map report is finite");
    }
}

void testConfigParsing() {
    const auto path = std::filesystem::temp_directory_path() / "cardillo_pj_termination_test.config";
    {
        std::ofstream configFile(path);
        configFile << "solver.criterion = NaTuRaL_MaP\n";
        configFile << "pj.final_check = Natural_Map\n";
        check((bool)configFile, "write temporary PJ config");
    }

    const config::Config cfg = config::ConfigReader::fromFile(path.string());
    std::error_code error;
    std::filesystem::remove(path, error);
    check(cfg.pj_criterion == config::PjTerminationCriterion::NaturalMap, "criterion config accepts solver alias and case-insensitive value");
    check(cfg.pj_final_check == config::PjFinalCheck::NaturalMap, "final-check config accepts case-insensitive value");
    check(config::Config{}.pj_criterion == config::PjTerminationCriterion::VelocityRms, "legacy RMS remains the default criterion");
    check(config::Config{}.pj_final_check == config::PjFinalCheck::None, "final-check defaults to disabled");
}

}  // namespace

int main() {
    testNaturalMapResidual();
    testBlockScaling();
    testTangentFrameInvariance();
    testNaturalMapContactSolve(false);
    testNaturalMapContactSolve(true);
    testConfigParsing();
    if (failures == 0) std::cout << "test_pj_termination: all checks passed\n";
    return failures == 0 ? 0 : 1;
}