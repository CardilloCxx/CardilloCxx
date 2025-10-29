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
        PhysicsSystem::Plane ground;
        sys.addObstacleBody(ground);

        const real_t mass = (real_t)1.0;
        const real_t radius = (real_t)0.1; // 10 cm radius
        
        // High-spin sphere dropped from height
        const Vector3r position1(0.0, 0.0, 1.0);
        const Vector3r angularVelocity1((real_t)0.0, (real_t)50.0, (real_t)0.0);
        sys.addRigidBodySphere(mass, position1, Quaternion4r::Identity(), Vector3r::Zero(), angularVelocity1, radius);

        // Lower-spin sphere
        const Vector3r position2(0.0, 1.0, 1.0);
        const Vector3r angularVelocity2((real_t)0.0, (real_t)10.0, (real_t)0.0);
        sys.addRigidBodySphere(mass, position2, Quaternion4r::Identity(), Vector3r::Zero(), angularVelocity2, radius);
    }
};
