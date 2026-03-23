#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

using namespace cardillo;

// Euler disk scene: approximate a spinning coin/cylinder on a plane using a
// short capsule (cylinder with rounded ends). The body is given a high
// initial spin and a slight tilt so it precesses on the plane.
class EulerDiskScene : public SceneBase {
public:
    EulerDiskScene() = default;
    ~EulerDiskScene() override = default;

    const char* sceneName() const override { return "euler_disk"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        // Ground plane
        physics::CubeShape groundShape(Vector3r(1, 1, 0.1));
        physics::RigidState groundState; groundState.position = Vector3r(0,0,-0.1); groundState.orientation = Quaternion4r::Identity();
        engine.addStaticBody(groundShape, groundState);

        // Cylinder approximation
        const real_t radius = (real_t)0.05;      // 5 cm radius
        const real_t halfLength = (real_t)0.01; // thin disk (2 cm thick)
        const real_t mass = (real_t)0.25;

        // Slight tilt so the disk contacts the plane at an edge and can precess
        const real_t tilt = (real_t) M_PI_2 + 0.4; // radians
        Quaternion4r q = Quaternion4r(Eigen::AngleAxis<real_t>(tilt, Vector3r::UnitX()));

        physics::RigidState state; 
        state.position = Vector3r(0.0, 0.0, radius + (real_t)0.001);
        state.orientation = q;
        state.angularVelocity = q.toRotationMatrix().transpose() * Vector3r((real_t)0.00, (real_t)0.0, (real_t)50.0);
        physics::RigidProps props(mass);
        entt::entity e = engine.addRigidBody(physics::CylinderShape(radius, halfLength), state, props);
        engine.track(e, "euler_disk_cylinder");
    }
};
