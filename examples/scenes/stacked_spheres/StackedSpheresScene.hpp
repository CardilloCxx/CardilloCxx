#pragma once
#include "../SceneBase.hpp"
#include <vector>
using namespace cardillo::physics;
class StackedSpheresScene : public SceneBase {
public:
    const char* sceneName() const override { return "stacked_spheres"; }
    StackedSpheresScene() = default;
    ~StackedSpheresScene() override = default;
    void populate(PhysicsEngine& engine) override {
        // Tube dimensions
        size_t nSpheres = 20;
        double sphereRadius = 0.05;
        double sphereSpacing = 0.0;

        bool disturb = false;

        double wallThickness = sphereRadius * 0.25;
        double tubeWidth = sphereRadius * 3.0;
        double tubeDepth = tubeWidth; 
        double tubeHeight = nSpheres * (2 * sphereRadius + sphereSpacing) + 0.5;

        double normalMass = 1.0;

        // Create 4 static cubes as tube walls
        auto w1 = engine.addStaticBody(CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}), RigidState{{-tubeWidth / 2.0 - wallThickness / 2.0, 0, tubeHeight / 2.0}});
        auto w2 = engine.addStaticBody(CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}), RigidState{{tubeWidth / 2.0 + wallThickness / 2.0, 0, tubeHeight / 2.0}});
        auto w3 = engine.addStaticBody(CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}), RigidState{{0, -tubeDepth / 2.0 - wallThickness / 2.0, tubeHeight / 2.0}});
        auto w4 = engine.addStaticBody(CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}), RigidState{{0, tubeDepth / 2.0 + wallThickness / 2.0, tubeHeight / 2.0}});

        // Place bottom cube (replaces bottom sphere)
        double bottomCubeHeight = 6 * sphereRadius;
        bottomCube = engine.addRigidBody(CubeShape({tubeWidth / 2.0, tubeDepth / 2.0, bottomCubeHeight / 2.0}), RigidState{{0.0, 0.0, bottomCubeHeight / 2.0}}, RigidProps(bottomMass));

        engine.disableCollisionBetween(bottomCube, w1);
        engine.disableCollisionBetween(bottomCube, w2);
        engine.disableCollisionBetween(bottomCube, w3);
        engine.disableCollisionBetween(bottomCube, w4);

        // cube_constraint = engine.addRigidConstraint(bottomCube);

        // Vertical velocity amplitude * sin(2 pi f t) and a constant spin of 5 rad/s about the body
        // x-axis, prescribed through the pose.
        const Vector3r p0(0.0, 0.0, bottomCubeHeight / 2.0);
        const real_t w = 2 * M_PI * frequency;
        engine.makeStatic(bottomCube);
        engine.addTrajectory(bottomCube, [=, this](real_t t) -> TrajectoryPose {
            const real_t z = amplitude / w * (1 - std::cos(w * t));
            return TrajectoryPose{p0 + Vector3r(0, 0, z), Quaternion4r(Eigen::AngleAxis<real_t>(5.0 * t, Vector3r::UnitX()))};
        });

        // Stack spheres above the bottom cube
        double z = bottomCubeHeight + sphereRadius;
        for (size_t i = 1; i < nSpheres; ++i) {
            engine.addRigidBody(SphereShape(sphereRadius), RigidState{{0, i * disturb * 1e-12, z}}, RigidProps(normalMass));
            z += 2 * sphereRadius + sphereSpacing;
        }
    }

    void updateScene(PhysicsEngine& /*engine*/, real_t /*t*/) override {}

private:
 double amplitude = 0.01;
 double frequency = 10.0;
 double bottomMass = 1.0;
 entt::entity bottomCube{entt::null};
 index_t cube_constraint{-1};
};
