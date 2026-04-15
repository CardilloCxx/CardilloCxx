#pragma once
#include "../SceneBase.hpp"
#include <vector>

class StackedSpheresScene : public SceneBase {
public:
    const char* sceneName() const override { return "stacked_spheres"; }
    StackedSpheresScene() = default;
    ~StackedSpheresScene() override = default;
    void populate(cardillo::physics::PhysicsEngine& engine) override {
        // Tube dimensions
        size_t nSpheres = 20;
        double sphereRadius = 0.05;
        double sphereSpacing = 0.0;

        bool disturb = true;

        double wallThickness = sphereRadius * 0.25;
        double tubeWidth = sphereRadius * 3.0;
        double tubeDepth = tubeWidth;
        double tubeHeight = nSpheres * (2 * sphereRadius + sphereSpacing) + 0.5;

        double normalMass = 1.0;

        // Create 4 static cubes as tube walls
        auto w1 = engine.addStaticBody(cardillo::physics::CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}),
                                       cardillo::physics::RigidState{{-tubeWidth / 2.0 - wallThickness / 2.0, 0, tubeHeight / 2.0}});
        auto w2 = engine.addStaticBody(cardillo::physics::CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}),
                                       cardillo::physics::RigidState{{tubeWidth / 2.0 + wallThickness / 2.0, 0, tubeHeight / 2.0}});
        auto w3 = engine.addStaticBody(cardillo::physics::CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}),
                                       cardillo::physics::RigidState{{0, -tubeDepth / 2.0 - wallThickness / 2.0, tubeHeight / 2.0}});
        auto w4 = engine.addStaticBody(cardillo::physics::CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}),
                                       cardillo::physics::RigidState{{0, tubeDepth / 2.0 + wallThickness / 2.0, tubeHeight / 2.0}});

        // Place bottom cube (replaces bottom sphere)
        double bottomCubeHeight = 6 * sphereRadius;
        bottomCube = engine.addRigidBody(cardillo::physics::CubeShape({tubeWidth / 2.0, tubeDepth / 2.0, bottomCubeHeight / 2.0}), cardillo::physics::RigidState{{0, 0, bottomCubeHeight / 2.0}},
                                         cardillo::physics::RigidProps(bottomMass));

        engine.disableCollisionBetween(bottomCube, w1);
        engine.disableCollisionBetween(bottomCube, w2);
        engine.disableCollisionBetween(bottomCube, w3);
        engine.disableCollisionBetween(bottomCube, w4);
        cube_constraint = engine.addRigidConstraint(bottomCube);

        // Stack spheres above the bottom cube
        double z = bottomCubeHeight + sphereRadius;
        for (size_t i = 1; i < nSpheres; ++i) {
            auto ent = engine.addRigidBody(cardillo::physics::SphereShape(sphereRadius), cardillo::physics::RigidState{{0, i * disturb * 1e-12, z}}, cardillo::physics::RigidProps(normalMass));
            sphereIds.push_back(ent);
            z += 2 * sphereRadius + sphereSpacing;
        }
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        double v = amplitude * std::sin(2 * M_PI * frequency * t);
        engine.setConstraintLinearVelocity(cube_constraint, Vector3r(0, 0, v));
    }
private:
    double amplitude = 0.5;
    double frequency = 0.6;
    double bottomMass = 1.0;
    std::vector<entt::entity> sphereIds;
    entt::entity bottomCube{entt::null};
    index_t cube_constraint{-1};
};
