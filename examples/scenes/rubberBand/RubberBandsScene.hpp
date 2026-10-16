#pragma once

#include "../SceneBase.hpp"
#include "misc/spline.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <random>

using namespace cardillo;

class RubberBandsScene : public SceneBase {
public:
    const char* sceneName() const override { return "rubber_bands"; }
    RubberBandsScene() = default;
    ~RubberBandsScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        engine.setGravity(Vector3r(0, 0, -9.81));
        engine.addStaticBody(physics::CubeShape(Vector3r(15.0, 15.0, 0.1)),
                             physics::RigidState(Vector3r(0, 0, -0.1)));

        const size_t numBands   = 12;
        const size_t segments   = 80;                // beam segments per band (~3 mm each)
        const real_t R          = (real_t)0.03;      // centerline radius (6 cm diameter)
        const real_t thickness  = (real_t)0.001;     // 1 mm, in the plane of the ring (radial)
        const real_t width      = (real_t)0.005;     // 5 mm, along the ring normal (flat band)
        const real_t E          = (real_t)1e6;       // soft rubber
        const real_t nu         = (real_t)0.49;      // nearly incompressible
        const real_t rho        = (real_t)1100.0;    // kg/m^3

        const auto section = physics::BeamCrossSection::rectangular(thickness, width);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);

        physics::RigidProps props = physics::RigidProps::withDensity(rho);
        props.friction = (real_t)0.8;

        misc::CircleSpline ring(Vector3r::Zero(), R);

        const real_t maxTilt      = (real_t)(25.0 * M_PI / 180.0);
        const real_t scatter      = (real_t)0.02;   // lateral jitter [m]
        const real_t firstHeight  = (real_t)0.03;   // height of the lowest band's center
        const real_t layerSpacing = (real_t)0.04;

        std::mt19937 rng(1234);
        std::uniform_real_distribution<real_t> uJitter(-scatter, scatter);
        std::uniform_real_distribution<real_t> uTilt((real_t)0, maxTilt);
        std::uniform_real_distribution<real_t> uPhi((real_t)0, (real_t)(2.0 * M_PI));

        for (size_t i = 0; i < numBands; ++i) {
            const real_t phi  = uPhi(rng);
            const real_t tilt = uTilt(rng);
            const Vector3r tiltAxis((real_t)std::cos(phi), (real_t)std::sin(phi), (real_t)0);
            const Quaternion4r q(Eigen::AngleAxis<real_t>(tilt, tiltAxis));
            const Vector3r pos(uJitter(rng), uJitter(rng), firstHeight + (real_t)i * layerSpacing);

            physics::RigidState state(pos, Vector3r::Zero(), q);

            auto ends = engine.createBeam(ring, section, springs, state, props, segments, physics::BeamColliderMode::InterSegmentHull);
        }
    }

    void updateScene(physics::PhysicsEngine& /*engine*/, real_t /*t*/, real_t /*dt*/) override {}
};
