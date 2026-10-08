// Consistency test for TranslationRotationConstraint: the force directions W must be the
// Jacobian of the position error, i.e. d/dt g = W_A^T u_A + W_B^T u_B for every row, including
// joints whose bodies have non-aligned frames in the reference configuration and large rotations
// about the joint axis. Checked by a central finite difference along random body velocities.

#include <cmath>
#include <iostream>
#include <random>

#include <Eigen/Geometry>

#include "physics/api/physics_engine.hpp"
#include "physics/constraints/constraints.hpp"
#include "rigid_body/rigid_body.hpp"

using namespace cardillo;
using namespace cardillo::physics;
using cardillo::RigidBody::RigidState;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        ++g_failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

Quaternion4r randomQuaternion(std::mt19937& rng) {
    std::normal_distribution<real_t> n(0, 1);
    Quaternion4r q(n(rng), n(rng), n(rng), n(rng));
    q.normalize();
    return q;
}

Vector3r randomVector(std::mt19937& rng, real_t scale = 1) {
    std::normal_distribution<real_t> n(0, 1);
    return scale * Vector3r(n(rng), n(rng), n(rng));
}

struct Pose {
    Vector3r r;
    Quaternion4r q;
};

void setPose(entt::registry& reg, entt::entity e, const Pose& p) {
    reg.get<C_Position3>(e).value = p.r;
    reg.get<C_Orientation>(e).setValue(p.q);
    RigidBody::updateState(reg, e);
}

// Drift of a pose with inertial linear velocity v and body-fixed angular velocity w over time s.
Pose drift(const Pose& p, const Vector3r& v, const Vector3r& w, real_t s) {
    Pose out;
    out.r = p.r + s * v;
    const real_t a = w.norm() * s;
    const Quaternion4r dq = (a != 0) ? Quaternion4r(Eigen::AngleAxis<real_t>(a, w.normalized())) : Quaternion4r::Identity();
    out.q = (p.q * dq).normalized();
    return out;
}

// Returns the largest error |FD(g) - W^T u| over all six rows, relative to max(1, |W^T u|).
real_t jacobianError(entt::registry& reg, physics::ConstraintPattern& c, entt::entity a, entt::entity b, const Pose& pa, const Pose& pb, std::mt19937& rng) {
    const Vector3r va = randomVector(rng), wa = randomVector(rng), vb = randomVector(rng), wb = randomVector(rng);
    const real_t eps = 1e-6;

    setPose(reg, a, drift(pa, va, wa, eps));
    setPose(reg, b, drift(pb, vb, wb, eps));
    const VectorXr gp = c.getConstraint().positionError;
    setPose(reg, a, drift(pa, va, wa, -eps));
    setPose(reg, b, drift(pb, vb, wb, -eps));
    const VectorXr gm = c.getConstraint().positionError;

    setPose(reg, a, pa);
    setPose(reg, b, pb);
    const physics::ConstraintResult res = c.getConstraint();
    Vector6r ua, ub;
    ua << va, wa;
    ub << vb, wb;
    const VectorXr Wu = res.WgA.transpose() * ua + res.WgB.transpose() * ub;
    const VectorXr fd = (gp - gm) / (2 * eps);
    return (fd - Wu).cwiseAbs().maxCoeff() / std::max((real_t)1, Wu.cwiseAbs().maxCoeff());
}

}  // namespace

int main() {
    std::mt19937 rng(42);

    for (int trial = 0; trial < 20; ++trial) {
        config::Config cfg;
        physics::PhysicsEngine engine(cfg);
        auto& reg = engine.ecs();

        // Two rigid bodies with arbitrary, non-aligned reference orientations.
        const Pose a0{randomVector(rng), randomQuaternion(rng)};
        const Pose b0{a0.r + randomVector(rng, 0.5), randomQuaternion(rng)};
        RigidProps props;
        props.mass = 1.0;
        const auto A = engine.addRigidBody(CubeShape(Vector3r(0.1, 0.2, 0.3)), RigidState{a0.r, a0.q}, props);
        const auto B = engine.addRigidBody(CubeShape(Vector3r(0.3, 0.1, 0.2)), RigidState{b0.r, b0.q}, props);
        setPose(reg, A, a0);
        setPose(reg, B, b0);

        const Vector3r jointPos = 0.5 * (a0.r + b0.r) + randomVector(rng, 0.1);
        const Vector3r axis = randomVector(rng).normalized();
        const size_t idx = engine.addTranslationRotationConstraint(A, B, physics::JointFrame::fromAxis(jointPos, axis), Vector3r::Constant(INFINITY), Vector3r::Zero(),
                                                                   Vector3r::Constant(INFINITY), Vector3r::Zero());
        auto& c = *engine.world().constraintPatterns()[idx];

        // 1) Zero error in the reference configuration.
        const VectorXr g0 = c.getConstraint().positionError;
        check(g0.cwiseAbs().maxCoeff() < 1e-12, "position error vanishes in the reference configuration (trial " + std::to_string(trial) + ")");

        // 2) Jacobian consistency at a nearby configuration (small relative motion).
        const Pose a1 = drift(a0, randomVector(rng, 0.01), randomVector(rng, 0.01), 1.0);
        const Pose b1 = drift(b0, randomVector(rng, 0.01), randomVector(rng, 0.01), 1.0);
        const real_t errNear = jacobianError(reg, c, A, B, a1, b1, rng);
        check(errNear < 1e-6, "W consistent with d/dt g near reference (trial " + std::to_string(trial) + ", err " + std::to_string(errNear) + ")");

        // 3) Jacobian consistency far from the reference configuration (large relative rotation).
        const Pose a2 = drift(a0, randomVector(rng, 0.3), randomVector(rng, 1.0), 1.0);
        const Pose b2 = drift(b0, randomVector(rng, 0.3), randomVector(rng, 1.0), 1.0);
        const real_t errFar = jacobianError(reg, c, A, B, a2, b2, rng);
        check(errFar < 1e-6, "W consistent with d/dt g far from reference (trial " + std::to_string(trial) + ", err " + std::to_string(errFar) + ")");

        // 4) Hinge: rotating B about the joint axis by a large angle leaves the two locked
        //    rotational rows (y, z) exactly zero.
        const Vector3r axisB = b0.q.conjugate() * axis;  // joint axis in B's body frame
        const Pose b3{jointPos + Eigen::AngleAxis<real_t>(2.5, axis) * (b0.r - jointPos), (b0.q * Quaternion4r(Eigen::AngleAxis<real_t>(2.5, axisB))).normalized()};
        setPose(reg, A, a0);
        setPose(reg, B, b3);
        const VectorXr g3 = c.getConstraint().positionError;
        check(std::abs(g3(4)) < 1e-12 && std::abs(g3(5)) < 1e-12 && g3.head<3>().cwiseAbs().maxCoeff() < 1e-12,
              "rotation about the joint axis leaves translational and locked rotational rows at zero (trial " + std::to_string(trial) + ")");
        check(std::abs(g3(3)) > 1e-3, "rotation about the joint axis is visible in the axis row (trial " + std::to_string(trial) + ")");

        // 5) Hinge: the force directions of the two locked rotational rows stay orthonormal for
        //    every hinge angle (no rank loss of the joint).
        for (real_t angle : {0.5, 1.5707963267948966, 2.5, 3.141592653589793}) {
            const Pose b4{jointPos + Eigen::AngleAxis<real_t>(angle, axis) * (b0.r - jointPos), (b0.q * Quaternion4r(Eigen::AngleAxis<real_t>(angle, axisB))).normalized()};
            setPose(reg, B, b4);
            const physics::ConstraintResult r4 = c.getConstraint();
            // inertial directions of rows y and z from body A's (static-free) block: n = A_IK1 * W_A,rot
            const Matrix33r RA = RigidBody::getState(reg, A).rotation;
            const Vector3r ny = RA * r4.WgA.block<3, 1>(3, 4);
            const Vector3r nz = RA * r4.WgA.block<3, 1>(3, 5);
            check(std::abs(ny.norm() - 1) < 1e-12 && std::abs(nz.norm() - 1) < 1e-12 && std::abs(ny.dot(nz)) < 1e-12,
                  "locked hinge rows keep orthonormal directions at angle " + std::to_string(angle) + " (trial " + std::to_string(trial) + ")");
        }
        setPose(reg, B, b0);
    }

    if (g_failures == 0) std::cout << "test_joint_jacobian: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
