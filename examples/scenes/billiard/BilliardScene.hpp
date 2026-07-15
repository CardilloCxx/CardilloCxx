#pragma once

#include "../SceneBase.hpp"

#include <vector>
#include <random>

using namespace cardillo;

class BilliardScene : public SceneBase {
   public:
    const char* sceneName() const override { return "billiard"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        engine.setGravity(Vector3r(0, 0, -9.81));

        // Simple rectangular table top.
        const Vector3r tableHalfExtents((real_t)0.99, (real_t)0.495, (real_t)0.05);
        const Vector3r tableCenter((real_t)0.0, (real_t)0.0, (real_t)0.71);
        engine.addStaticBody(physics::CubeShape(tableHalfExtents), physics::RigidState(tableCenter));

        const real_t railThickness = (real_t)0.02;
        const real_t railZ = tableCenter.z() + tableHalfExtents.z() + (real_t)0.02;

        const Vector3r longRailExtents(tableHalfExtents.x(), railThickness, (real_t)0.05);
        engine.addStaticBody(physics::CubeShape(longRailExtents), physics::RigidState(Vector3r(0.0, tableHalfExtents.y() + railThickness, railZ)));
        engine.addStaticBody(physics::CubeShape(longRailExtents), physics::RigidState(Vector3r(0.0, -(tableHalfExtents.y() + railThickness), railZ)));

        const Vector3r shortRailExtents(railThickness, tableHalfExtents.y(), (real_t)0.05);
        engine.addStaticBody(physics::CubeShape(shortRailExtents), physics::RigidState(Vector3r(tableHalfExtents.x() + railThickness, 0.0, railZ)));
        engine.addStaticBody(physics::CubeShape(shortRailExtents), physics::RigidState(Vector3r(-(tableHalfExtents.x() + railThickness), 0.0, railZ)));

        const real_t ballRadius = (real_t)0.0285;
        const real_t ballSpacing = (real_t)2.0 * ballRadius;
        const real_t tableTopZ = tableCenter.z() + tableHalfExtents.z() + ballRadius;

        const Vector3r rackOrigin((real_t)0.0, (real_t)0.0, tableTopZ);
        const int rows = 5;
        const real_t rowSpacingY = ballSpacing * (real_t)0.8660254037844386; // sqrt(3)/2

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(-0.001f, 0.001f);

        for (int row = 0; row < rows; ++row) {
            const int count = row + 1;
            const real_t rowOffset = (real_t)(count - 1) * ballSpacing * (real_t)0.5;
            const real_t rowY = rackOrigin.y() + (real_t)row * rowSpacingY;

            for (int col = 0; col < count; ++col) {
                const real_t x = rackOrigin.x() + (real_t)col * ballSpacing - rowOffset;

                const Vector3r pos(rowY + dist(gen), x + dist(gen), tableTopZ);

                physics::RigidProps props = physics::RigidProps::withDensity((real_t)1700.0);
                props.friction = (real_t)0.15;
                props.restitution_normal = (real_t)0.9;
                props.restitution_tangential = (real_t)0.1;
                engine.addRigidBody(physics::SphereShape(ballRadius), physics::RigidState(pos), props);
            }
        }

        // Cue ball placed in front of the rack and given an initial push.
        physics::RigidState cueState;
        cueState.position = Vector3r((real_t)-0.20, (real_t)0.0 + dist(gen), tableTopZ);
        cueState.linearVelocity = Vector3r((real_t)9.0 + dist(gen), (real_t)0.0 + dist(gen), (real_t)0.0 + dist(gen));
        physics::RigidProps cueProps = physics::RigidProps::withDensity((real_t)1700.0);
        cueProps.friction = (real_t)0.15;
        cueProps.restitution_normal = (real_t)0.9;
        cueProps.restitution_tangential = (real_t)0.1;
        engine.addRigidBody(physics::SphereShape(ballRadius), cueState, cueProps);
    }
};
