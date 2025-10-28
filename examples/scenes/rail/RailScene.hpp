#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

// Rail scene: loads a rail mesh as a static obstacle, a ground plane, and a
// train wheel as a dynamic rigid mesh. The wheel receives an initial linear
// velocity along the rail and an angular velocity to simulate rolling.
class RailScene : public SceneBase {
public:
    RailScene() = default;
    ~RailScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane (thin cube)
        PhysicsSystem::Cube ground;
        const real_t groundHalfThickness = (real_t)0.01;
        const real_t groundHalfSize = (real_t)50.0;
        ground.center = Vector3r(0.0, 0.0, -groundHalfThickness);
        ground.halfExtents = Vector3r(groundHalfSize, groundHalfSize, groundHalfThickness);
        ground.q = Quaternion4r::Identity();
        sys.addObstacleBody(ground);

        // Rail mesh as a static obstacle (long rail)
        const std::string railMesh = "res/meshes/rail.obj";
        const Vector3r railPos(0.0, 0.0, 0.0);
        const Quaternion4r railOri = Quaternion4r::Identity();
        const Vector3r railScale(1.0, 1.0, 1.0);
        sys.addObstacleMesh(railPos, railOri, railMesh, railScale);

        // Wheel mesh as a dynamic rigid body
        const std::string wheelMesh = "res/meshes/train_wheel.obj";
        const Vector3r wheelPos(-0.075, 0.0, 0.575);
        const Quaternion4r wheelOri =  Quaternion4r::Identity();
        const Vector3r wheelScale(1.0, 1.0, 1.0);
        const real_t wheelMass = (real_t)5.0;
        const Vector3r linearVel(0.0, -1.0, 0.0);
        const real_t wheelRadius = (real_t)0.5;
        const real_t trackRadius = (real_t)10.0;
        const Vector3r angularVel(-linearVel.y() / wheelRadius, 0.0,  -linearVel.y() / trackRadius);

        sys.addRigidBodyMesh(wheelMass, wheelPos, wheelOri, linearVel, angularVel, wheelMesh, wheelScale);
    }
};
