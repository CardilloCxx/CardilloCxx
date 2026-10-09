// Consistency test for BeamConstraint (rod element) force directions W versus the time derivative
// of the strain measures g = (gamma, kappa). The mixed Petrov--Galerkin directions (Dai et al. 2026,
// Eq. (20)) are derived from the internal virtual work with linearly interpolated virtual rotations;
// they are not the Jacobian of g but coincide with it up to second order in the relative rotation
// phi of the two nodes.

#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Geometry>

#include "physics/api/physics_engine.hpp"
#include "physics/constraints/constraints.hpp"
#include "rigid_body/rigid_body.hpp"

using namespace cardillo;
using namespace cardillo::physics;

namespace {

struct Pose {
    Vector3r r;
    Quaternion4r q;
};

void setPose(entt::registry& reg, entt::entity e, const Pose& p) {
    reg.get<C_Position3>(e).value = p.r;
    reg.get<C_Orientation>(e).setValue(p.q);
    RigidBody::updateState(reg, e);
}

Pose drift(const Pose& p, const Vector3r& v, const Vector3r& w, real_t s) {
    const real_t a = w.norm() * s;
    const Quaternion4r dq = (a != 0) ? Quaternion4r(Eigen::AngleAxis<real_t>(a, w.normalized())) : Quaternion4r::Identity();
    return {p.r + s * v, (p.q * dq).normalized()};
}

Vector3r randomVector(std::mt19937& rng, real_t scale = 1) {
    std::normal_distribution<real_t> n(0, 1);
    return scale * Vector3r(n(rng), n(rng), n(rng));
}

// Returns |FD(g) - W^T u| / max(1, |FD(g)|) for a random state; phi is the relative rotation angle.
real_t jacobianDeviation(real_t spread, std::mt19937& rng, real_t& phi) {
    config::Config cfg;
    PhysicsEngine engine(cfg);
    auto& reg = engine.ecs();

    std::normal_distribution<real_t> n(0, 1);
    const Quaternion4r qa = Quaternion4r(n(rng), n(rng), n(rng), n(rng)).normalized();
    const Quaternion4r qb = (qa * Quaternion4r(Eigen::AngleAxis<real_t>(spread, randomVector(rng).normalized()))).normalized();
    const Pose a0{randomVector(rng), qa};
    const Pose b0{a0.r + qa * Vector3r(0.1, 0.01, -0.02), qb};

    RigidProps props;
    props.mass = 1.0;
    const auto A = engine.addRigidBody(CubeShape(Vector3r(0.05, 0.01, 0.01)), RigidState{a0.r, a0.q}, props);
    const auto B = engine.addRigidBody(CubeShape(Vector3r(0.05, 0.01, 0.01)), RigidState{b0.r, b0.q}, props);
    setPose(reg, A, a0);
    setPose(reg, B, b0);
    BeamSpringParams springs = BeamSpringParams::fromMaterial(2e11, 0.3);
    const size_t idx = engine.addBeamConstraint(A, B, springs, BeamCrossSection(0.01, 0.01, BeamBodyType::Capsule));
    auto& c = *engine.world().constraintPatterns()[idx];

    // Evaluate at a configuration away from the reference (as during a simulation).
    const Pose a1 = drift(a0, randomVector(rng, 0.01), randomVector(rng, 0.2 * spread), 1.0);
    const Pose b1 = drift(b0, randomVector(rng, 0.01), randomVector(rng, 0.2 * spread), 1.0);
    phi = 2 * std::acos(std::min((real_t)1, std::abs(a1.q.dot(b1.q))));
    const Vector3r va = randomVector(rng), wa = randomVector(rng), vb = randomVector(rng), wb = randomVector(rng);
    const real_t eps = 1e-6;
    setPose(reg, A, drift(a1, va, wa, eps));
    setPose(reg, B, drift(b1, vb, wb, eps));
    const VectorXr gp = c.getConstraint().positionError;
    setPose(reg, A, drift(a1, va, wa, -eps));
    setPose(reg, B, drift(b1, vb, wb, -eps));
    const VectorXr gm = c.getConstraint().positionError;
    setPose(reg, A, a1);
    setPose(reg, B, b1);
    const ConstraintResult res = c.getConstraint();
    Vector6r ua, ub;
    ua << va, wa;
    ub << vb, wb;
    const VectorXr Wu = res.WgA.transpose() * ua + res.WgB.transpose() * ub;
    const VectorXr fd = (gp - gm) / (2 * eps);
    return (fd - Wu).cwiseAbs().maxCoeff() / std::max((real_t)1, fd.cwiseAbs().maxCoeff());
}

}  // namespace

int main() {
    std::mt19937 rng(7);
    int failures = 0;
    real_t worstRatio = 0;
    for (int trial = 0; trial < 20; ++trial) {
        const real_t spread = 0.05 + 0.45 * (trial % 10) / 9.0;
        real_t phi = 0;
        const real_t err = jacobianDeviation(spread, rng, phi);
        worstRatio = std::max(worstRatio, err / (phi * phi));
        if (err > phi * phi) {
            ++failures;
            std::cerr << "FAILED: Petrov-Galerkin directions not consistent to first order, trial " << trial << " (phi " << phi << "): deviation " << err << std::endl;
        }
    }
    std::cout << "test_beam_jacobian: Petrov-Galerkin worst deviation / phi^2 " << worstRatio << std::endl;
    if (failures == 0) std::cout << "test_beam_jacobian: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
