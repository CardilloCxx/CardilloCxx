#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

// Rotating ball scene: a single sphere dropped onto a floor with high
// initial angular velocity (w = 50 rad/s) to exercise spinning/friction.
class RotatingBallScene : public SceneBase {
public:
    RotatingBallScene() = default;
    ~RotatingBallScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane for visual + collision
    PhysicsSystem::PlaneShape groundShape{Vector3r(0,0,1), Vector3r(0,1,0), (real_t)5.0, (real_t)5.0};
    PhysicsSystem::RigidState groundState; groundState.position = Vector3r::Zero(); groundState.orientation = Quaternion4r::Identity();
    sys.addStaticBody(groundShape, groundState);

        const real_t mass = (real_t)1.0;
        const real_t radius = (real_t)0.1; // 10 cm radius
        
        // High-spin sphere dropped from height
        const Vector3r position1(0.0, 0.0, 1.0);
        const Vector3r angularVelocity1((real_t)0.0, (real_t)50.0, (real_t)0.0);
        {
            sys.addRigidBody(PhysicsSystem::SphereShape(radius), PhysicsSystem::RigidState(position1, Vector3r::Zero(), angularVelocity1), PhysicsSystem::RigidProps(mass));
        }

        // Lower-spin sphere
        const Vector3r position2(0.0, 1.0, 1.0);
        const Vector3r angularVelocity2((real_t)0.0, (real_t)10.0, (real_t)0.0);
        {
            sys.addRigidBody(PhysicsSystem::SphereShape(radius), PhysicsSystem::RigidState(position2, Vector3r::Zero(), angularVelocity2), PhysicsSystem::RigidProps(mass));
        }
    }
};
