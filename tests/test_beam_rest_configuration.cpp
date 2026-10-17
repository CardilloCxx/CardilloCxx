#include <cmath>
#include <iostream>

#include "physics/api/physics_engine.hpp"

using namespace cardillo;
using namespace cardillo::physics;

namespace {

int failures = 0;

void check(bool ok, const char* message) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

VectorXr beamError(const Vector3r& endpoint, const BeamSpringParams& springs) {
    config::Config config;
    PhysicsEngine engine(config);
    RigidProps props;
    props.mass = (real_t)1;
    const auto a = engine.addRigidBody(CubeShape(Vector3r::Constant((real_t)0.1)), RigidState(Vector3r::Zero(), Quaternion4r::Identity()), props);
    const auto b = engine.addRigidBody(CubeShape(Vector3r::Constant((real_t)0.1)), RigidState(endpoint, Quaternion4r::Identity()), props);
    const auto index = engine.addBeamConstraint(a, b, springs, BeamCrossSection::rectangular((real_t)0.1, (real_t)0.1));
    return engine.world().constraintPatterns()[index]->getConstraint().positionError;
}

}  // namespace

int main() {
    const auto initialError = beamError(Vector3r((real_t)1, (real_t)1, (real_t)0), BeamSpringParams{});
    check(initialError.norm() < (real_t)1e-12, "default rest state is relaxed at the initial pose");

    BeamSpringParams straight;
    straight.relaxedInStraightConfiguration();
    const auto straightError = beamError(Vector3r((real_t)1, (real_t)1, (real_t)0), straight);
    const Vector3r expectedStraightError((real_t)1 - std::sqrt((real_t)2), (real_t)1, (real_t)0);
    check((straightError.head<3>() - expectedStraightError).norm() < (real_t)1e-12, "straight rest shape removes initial shear");
    check(straightError.tail<3>().norm() < (real_t)1e-12, "straight rest shape has zero curvature");

    BeamSpringParams pretension;
    pretension.withAxialPretension((real_t)0.04);
    const auto lowStrainError = beamError(Vector3r((real_t)1, (real_t)0, (real_t)0), pretension);
    check(std::abs(lowStrainError.x() - ((real_t)1 - (real_t)1 / (real_t)1.04)) < (real_t)1e-12, "small pretension uses engineering strain");

    pretension.withAxialPretension((real_t)1.2);
    const auto highStrainError = beamError(Vector3r((real_t)1, (real_t)0, (real_t)0), pretension);
    check(std::abs(highStrainError.x() - ((real_t)1 - (real_t)1 / (real_t)2.2)) < (real_t)1e-12, "pretension strains greater than 1 are supported");

    pretension.withAxialPretension((real_t)-0.2);
    const auto compressionError = beamError(Vector3r((real_t)1, (real_t)0, (real_t)0), pretension);
    check(std::abs(compressionError.x() - ((real_t)1 - (real_t)1 / (real_t)0.8)) < (real_t)1e-12, "negative strain represents a compressed initial state");

    if (failures == 0) std::cout << "test_beam_rest_configuration: all checks passed\n";
    return failures == 0 ? 0 : 1;
}