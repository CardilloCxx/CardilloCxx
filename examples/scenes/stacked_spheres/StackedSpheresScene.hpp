#pragma once
#include "../SceneBase.hpp"
#include <vector>

class StackedSpheresScene : public SceneBase {
public:
    const char* sceneName() const override { return "stacked_spheres"; }
    StackedSpheresScene() = default;
    ~StackedSpheresScene() override = default;
    void populate(cardillo::PhysicsSystem& sys) override {
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
        cardillo::physics::BodyFactory::addStaticBody(sys, cardillo::PhysicsSystem::CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}),
            cardillo::PhysicsSystem::RigidState{{-tubeWidth/2.0-wallThickness/2.0, 0, tubeHeight/2.0}});
        cardillo::physics::BodyFactory::addStaticBody(sys, cardillo::PhysicsSystem::CubeShape({wallThickness / 2, tubeDepth / 2, tubeHeight / 2}),
            cardillo::PhysicsSystem::RigidState{{tubeWidth/2.0+wallThickness/2.0, 0, tubeHeight/2.0}});
        cardillo::physics::BodyFactory::addStaticBody(sys, cardillo::PhysicsSystem::CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}),
            cardillo::PhysicsSystem::RigidState{{0, -tubeDepth/2.0-wallThickness/2.0, tubeHeight/2.0}});
        cardillo::physics::BodyFactory::addStaticBody(sys, cardillo::PhysicsSystem::CubeShape({tubeWidth / 2 + wallThickness, wallThickness / 2, tubeHeight / 2}),
            cardillo::PhysicsSystem::RigidState{{0, tubeDepth/2.0+wallThickness/2.0, tubeHeight/2.0}});

        // Place bottom cube (replaces bottom sphere)
        double bottomCubeHeight = 2 * sphereRadius;
        auto bottomCube = cardillo::physics::BodyFactory::addRigidBody(sys, 
            cardillo::PhysicsSystem::CubeShape({tubeWidth/2.0, tubeDepth/2.0, bottomCubeHeight/2.0}),
            cardillo::PhysicsSystem::RigidState{{0, 0, bottomCubeHeight/2.0}},
            cardillo::PhysicsSystem::RigidProps(bottomMass)
        );
        bottomSphere = bottomCube;

        // Stack spheres above the bottom cube
        double z = bottomCubeHeight + sphereRadius;
        for (size_t i = 1; i < nSpheres; ++i) {
            auto ent = cardillo::physics::BodyFactory::addRigidBody(sys, 
                cardillo::PhysicsSystem::SphereShape(sphereRadius),
                cardillo::PhysicsSystem::RigidState{{0, i * disturb * 1e-12, z}},
                cardillo::PhysicsSystem::RigidProps(normalMass)
            );
            sphereIds.push_back(ent);
            z += 2 * sphereRadius + sphereSpacing;
        }
    }
    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        double v = amplitude * std::sin(2 * M_PI * frequency * t);
        sys.applyForce(bottomSphere, {0, 0, 9.81 * bottomMass}, {0, 0, 0});
        sys.setLinearVelocity(bottomSphere, {0, 0, v});
        sys.setAngularVelocity(bottomSphere, {0, 0, 0});
    }
private:
    double amplitude = 0.5;
    double frequency = 0.3;
    double bottomMass = 1e14;
    std::vector<entt::entity> sphereIds;
    entt::entity bottomSphere{entt::null};
};
