#pragma once

#include <cmath>
#include "../SceneBase.hpp"
#include "misc/spline.hpp"

using namespace cardillo;

class StringTestScene : public SceneBase {
   public:
    const char* sceneName() const override { return "StringTest"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        engine.setGravity(Vector3r(0, 0, -9.81));

        const Vector3r topCubeHalfExtents((real_t)0.22, (real_t)0.22, (real_t)0.03);
        const Vector3r topCubeCenter((real_t)0.0, (real_t)0.0, (real_t)0.45);
        m_topCube = engine.addStaticBody(physics::CubeShape(topCubeHalfExtents), physics::RigidState(topCubeCenter));

        const Vector3r floorHalfExtents((real_t)1.2, (real_t)1.2, (real_t)0.05);
        engine.addStaticBody(physics::CubeShape(floorHalfExtents), physics::RigidState(Vector3r(0, 0, -floorHalfExtents.z())));

        const int nx = 9;
        const int ny = 7;
        const real_t spacing = (real_t)0.045;
        const real_t zTop = topCubeCenter.z() - topCubeHalfExtents.z() - (real_t)0.002;
        const real_t stringLength = (real_t)0.26;
        const real_t zBottom = zTop - stringLength;
        const size_t segments = 16;

        physics::RigidProps stringProps = physics::RigidProps::withDensity((real_t)850.0);
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());

        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                const real_t x = ((real_t)ix - (real_t)(nx - 1) * (real_t)0.5) * spacing;
                const real_t y = ((real_t)iy - (real_t)(ny - 1) * (real_t)0.5) * spacing;

                const real_t u = (nx > 1) ? (real_t)ix / (real_t)(nx - 1) : (real_t)0;
                const real_t v = (ny > 1) ? (real_t)iy / (real_t)(ny - 1) : (real_t)0;

                const real_t d = (real_t)0.002 + (real_t)0.004 * v;
                const real_t E = (real_t)4.0e4 + (real_t)1.2e5 * u;
                const real_t nu = (real_t)0.35;
                const real_t damping = (real_t)0.04 + (real_t)0.26 * ((u + v) * (real_t)0.5);

                physics::BeamCrossSection section(d, d, physics::BeamBodyType::Capsule);
                auto springs = physics::BeamSpringParams::fromMaterial(E, nu, (real_t)50, (real_t)1, (real_t)1, (real_t)1, (real_t)1, damping);

                LinearSpline spline(Vector3r(x, y, zTop), Vector3r(x, y, zBottom));
                auto ends = engine.createBeam(spline, section, springs, stateDefaults, stringProps, segments);
                engine.makeStatic(ends.first);
                Vector3r pos = engine.getPosition(ends.second).head(3);
                auto weight = engine.addRigidBody(physics::SphereShape(spacing * 0.25), physics::RigidState(pos), physics::RigidProps::withDensity((real_t)200.0));
                engine.addRigidConstraint(ends.second, weight);

                engine.disableCollisionBetween(ends.first, m_topCube);
            }
        }

        const Vector3r combHalfExtents((real_t)0.014, (real_t)0.18, (real_t)0.024);
        const Vector3r combStart((real_t)-0.20, (real_t)0.0, zBottom + (real_t)0.08);

        m_combCube = engine.addRigidBody(physics::CubeShape(combHalfExtents), physics::RigidState(combStart, Vector3r::Zero(), Quaternion4r::Identity()), physics::RigidProps((real_t)45.0));

        engine.addTrajectory(m_combCube, std::nullopt, [](real_t t) {
            TrajectoryTwist twist;
            twist.first = Vector3r(0.5 * std::sin(t), 0, 0);
            twist.second = Vector3r::Zero();
            return twist;
        });
    }

   private:
    entt::entity m_topCube{entt::null};
    entt::entity m_combCube{entt::null};
};
