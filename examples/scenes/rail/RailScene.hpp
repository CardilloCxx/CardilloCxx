#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

// Rail scene: loads a rail mesh as a static obstacle, a ground plane, and a
// train wheel as a dynamic rigid mesh. The wheel receives an initial linear
// velocity along the rail and an angular velocity to simulate rolling.
class RailScene : public SceneBase {
public:
    const char* sceneName() const override { return "rail"; }
    RailScene() = default;
    ~RailScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane (thin cube)
    const real_t groundHalfThickness = (real_t)0.01;
    const real_t groundHalfSize = (real_t)50.0;
    PhysicsSystem::CubeShape groundShape{Vector3r(groundHalfSize, groundHalfSize, groundHalfThickness)};
    PhysicsSystem::RigidState groundState; groundState.position = Vector3r(0.0, 0.0, -groundHalfThickness); groundState.orientation = Quaternion4r::Identity();
    sys.addStaticBody(groundShape, groundState);

        // Rail mesh as a static obstacle (long rail)
    sys.addStaticBody(PhysicsSystem::MeshShape("res/meshes/rail.obj"), PhysicsSystem::RigidState(Vector3r(0,0,0))); // static

        // Wheel mesh as a dynamic rigid body
        const Vector3r linearVel(0.0, -1.0, 0.0);
        const real_t wheelRadius = (real_t)0.5;
        const real_t trackRadius = (real_t)10.0;
        const Vector3r angularVel(-linearVel.y() / wheelRadius, 0.0,  -linearVel.y() / trackRadius);
    PhysicsSystem::RigidState s(Vector3r(-0.075, 0.0, 0.575), linearVel, Quaternion4r::Identity(), angularVel);
    sys.addRigidBody(PhysicsSystem::MeshShape("res/meshes/train_wheel.obj"), s, PhysicsSystem::RigidProps((real_t)5.0));
    }
};
