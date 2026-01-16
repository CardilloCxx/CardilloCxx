#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace cardillo;

// Unilateral rocking-rod scene:
// Two static supports at x=±a and a tilted rod initially in point contact at +a.
class UnilateralScene : public SceneBase {
public:
    const char* sceneName() const override { return "unilateral"; }

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        const auto formatAForName = [](real_t aVal) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2) << (double)aVal;
            std::string s = ss.str();
            for (char& c : s) if (c == '.') c = 'p';
            return s;
        };

        const real_t l = (real_t)1.0;
        const real_t m = (real_t)1.0;

        const real_t phiDeg = (real_t)30.0;
        const real_t phi = phiDeg * (real_t)M_PI / (real_t)180.0; // rotation about +Y

        const std::vector<real_t> aValues = { (real_t)0.1, (real_t)0.2, (real_t)0.3 };
        const real_t ySpacing = (real_t)0.50; // separate copies to avoid cross-interactions

        // Keep peg/support geometry constant across variants (do not scale with a).
        // Baseline chosen to match the original a=0.1 setup.
        const real_t supportRadius0 = (real_t)0.005;
        const real_t supportWidth0  = (real_t)0.10; // total capsule length along Y
        const real_t supportZ0      = (real_t)0.02; // height of support center
        const real_t supportHalfLength0 = std::max<real_t>(
            (real_t)0.0,
            (supportWidth0 - (real_t)2.0 * supportRadius0) * (real_t)0.5
        );

        // Shared orientations
        const Quaternion4r qSupport(Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitY()));
        const Quaternion4r qCapZX(Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitX()));
        const Quaternion4r qRodTilt(Eigen::AngleAxis<real_t>(phi, Vector3r::UnitY()));
        const Quaternion4r qRod = qRodTilt * qCapZX;

        for (size_t i = 0; i < aValues.size(); ++i) {
            const real_t a = aValues[i];
            const Vector3r offset((real_t)0.0, (real_t)i * ySpacing, (real_t)0.0);

            // 1) Static supports: capsule-shaped pegs (constant size)
            const real_t supportRadius = supportRadius0;
            const real_t supportHalfLength = supportHalfLength0;
            const real_t H = supportZ0;
            sys.addStaticBody(
                PhysicsSystem::CapsuleShape(supportRadius, supportHalfLength),
                PhysicsSystem::RigidState(Vector3r(-a, (real_t)0.0, H) + offset, qSupport)
            );
            sys.addStaticBody(
                PhysicsSystem::CapsuleShape(supportRadius, supportHalfLength),
                PhysicsSystem::RigidState(Vector3r(+a, (real_t)0.0, H) + offset, qSupport)
            );

            // 2) Dynamic rocking rod: capsule-shaped, mass m.
            // Keep rod thickness constant across variants (do not scale with a).
            // Matches the a=0.1 baseline thickness: 0.05 * 0.1 = 0.005.
            const real_t rodRadius = (real_t)0.005;
            const real_t rodHalfLength = l * (real_t)0.5;

            // 3) Initial CoM position.
            const real_t x0 = (real_t)0.0;
            const real_t z_contact = H + supportRadius * (real_t)2.0; // top point of the +a support capsule
            const real_t z0 = z_contact + a * std::tan(phi);
            const Vector3r xcm0(x0, (real_t)0.0, z0);

            PhysicsSystem::RigidState rodState;
            rodState.position = xcm0 + offset;
            rodState.orientation = qRod;
            rodState.linearVelocity = Vector3r::Zero();
            rodState.angularVelocity = Vector3r::Zero();

            PhysicsSystem::RigidProps rodProps;
            rodProps.mass = m;

            auto rod = sys.addRigidBody(PhysicsSystem::CapsuleShape(rodRadius, rodHalfLength), rodState, rodProps);
            const std::string aTag = formatAForName(a);
            sys.track(rod, "unilateral_rod_a" + aTag);

            std::cout << "[Unilateral] a=" << a << ", m=" << m
                      << ", phi=" << phiDeg << " deg, track=unilateral_rod_a" << aTag
                      << ", offset=" << offset.transpose() << std::endl;
        }
    }
};
