#pragma once

#include <Eigen/Geometry>
#include <cmath>
#include "../SceneBase.hpp"

// Tube scene: approximate the inner wall of a cylinder with N static box
// segments and shoot a single sphere through one open end.
class TubeScene : public SceneBase {
   public:
    const char* sceneName() const override { return "tube"; }
    TubeScene() = default;
    ~TubeScene() override = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        engine.setGravity(Vector3r::Zero());

        const int nWallBoxes = 128;
        const real_t tubeLength = (real_t)0.3;
        const real_t innerRadius = (real_t)0.10;
        const real_t wallThickness = (real_t)0.01;
        const real_t halfLength = (real_t)0.5 * tubeLength;
        const real_t halfRadial = (real_t)0.5 * wallThickness;
        const real_t halfTangential = (real_t)1.01 * (innerRadius + wallThickness) * std::tan(M_PI / (real_t)nWallBoxes);

        for (int i = 0; i < nWallBoxes; ++i) {
            const real_t theta = ((real_t)2.0 * M_PI * (real_t)i) / (real_t)nWallBoxes;
            const real_t cy = (innerRadius + halfRadial) * std::cos(theta);
            const real_t cz = (innerRadius + halfRadial) * std::sin(theta);

            const Vector3r center((real_t)0.0, cy, cz);
            const Quaternion4r orientation(Eigen::AngleAxis<real_t>(theta, Vector3r::UnitX()));
            const Vector3r halfExtents(halfLength, halfRadial, halfTangential);

            engine.addStaticBody(physics::CubeShape{halfExtents}, physics::RigidState(center, orientation));
        }

        const real_t sphereRadius = (real_t)0.01;
        const Vector3r sphereStart(-(halfLength - sphereRadius), innerRadius - sphereRadius * 1.5, (real_t)0.0);
        const Vector3r sphereVelocity((real_t)2.5, (real_t)0.8, (real_t)2.55);

        entt::entity sphere = engine.addRigidBody(physics::SphereShape(sphereRadius), physics::RigidState(sphereStart, sphereVelocity), physics::RigidProps::withDensity((real_t)1000.0));

        engine.track(sphere, "tube_sphere");
    }
};
