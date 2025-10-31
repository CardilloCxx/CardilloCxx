#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include "physics/force_interaction.hpp"

using namespace cardillo;

class SpringTestScene : public SceneBase {
public:
    SpringTestScene() = default;
    ~SpringTestScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // ground (static obstacle cube)
        PhysicsSystem::Cube ground;
        ground.center = Vector3r(0.0, 0.0, -0.5);
        ground.halfExtents = Vector3r(5.0, 5.0, 0.5);
        ground.q = Quaternion4r::Identity();
        // addObstacleBody returns an integer id; convert to entt::entity to reference in the spring
        index_t ground_id = sys.addObstacleBody(ground);
        entt::entity eGround = entt::entity(static_cast<uint32_t>(ground_id));

        // single dynamic cube above the ground
        PhysicsSystem::Cube c; c.center = Vector3r(0.0, 0.0, 0.6); c.halfExtents = Vector3r(0.1,0.1,0.1); c.q = Quaternion4r::Identity();
        entt::entity eBox = sys.addRigidBody((real_t)1.0, c.center, c.q, Vector3r(0.0, 0.0, -20.0), Vector3r::Zero(), c);

        // Create a spring between the dynamic box (eBox) and the ground (eGround).
        cardillo::physics::SpringConstraint spring(sys.ecs(), eBox, eGround);
        spring.addTranslationalSpring(Vector3r::UnitZ(), (real_t)1000.0, (real_t)0.0);
        sys.addConstraint(spring);
    }
};
