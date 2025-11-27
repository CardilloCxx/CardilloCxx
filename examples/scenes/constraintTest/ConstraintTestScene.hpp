#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class ConstraintTestScene : public SceneBase {
public:
    const char* sceneName() const override { return "constraint_test"; }
    ConstraintTestScene() = default;
    ~ConstraintTestScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Floor (static)
        // sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(15.0, 15.0, 0.1)), PhysicsSystem::RigidState(Vector3r(0,0,-0.1)));
        const real_t cubeHalfExt = 0.5;

        // Cube with fixed constraint
        m_cube = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(cubeHalfExt, cubeHalfExt, cubeHalfExt)),
            PhysicsSystem::RigidState(Vector3r(0,0,2)),
            PhysicsSystem::RigidProps(1e10)
        );

        // small cube allowed to move in a plane
        auto cube2 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, 2 + cubeHalfExt + 0.05)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, cube2);
        sys.disableCollisionBetween(m_cube, cube2);

        // small cube allowed to move in a plane
        auto cube3 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(-cubeHalfExt - 0.05, cubeHalfExt + 0.05, 2 + cubeHalfExt + 0.05)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        auto jp = physics::JointProperties( sys.ecs(), m_cube, cube3);
        Vector3r K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(2) = 1000;
        Vector3r D_trans = Vector3r::Zero();
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube3, jp, K_trans, D_trans);
        sys.disableCollisionBetween(m_cube, cube3);

        auto cube4 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(-cubeHalfExt - 0.05, -cubeHalfExt - 0.05, 2 + cubeHalfExt + 0.05)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        jp = physics::JointProperties( sys.ecs(), m_cube, cube4);
        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(1) = 1000;
        D_trans = Vector3r::Zero();
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube4, jp, K_trans, D_trans);
        sys.disableCollisionBetween(m_cube, cube4);

        auto cube5 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, -cubeHalfExt - 0.05, 2 + cubeHalfExt + 0.05)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        jp = physics::JointProperties( sys.ecs(), m_cube, cube5);
        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(0) = 1000;
        D_trans = Vector3r::Zero();
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube5, jp, K_trans, D_trans);
        sys.disableCollisionBetween(m_cube, cube5);

         auto cube6 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, 2 - cubeHalfExt - 0.05)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        jp = physics::JointProperties( sys.ecs(), m_cube, cube6, std::nullopt, Vector3r(0, 0, 2), Vector3r(cubeHalfExt, cubeHalfExt, 2-cubeHalfExt), Matrix33r::Identity());
        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        D_trans = Vector3r::Zero();
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube6, jp, K_trans, D_trans);
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        // sys.setAngularVelocity(m_cube, Vector3r(0.1432, 0.0234, 0.0432)); 
        sys.setAngularVelocity(m_cube, Vector3r(0,0.4,0));
        // sys.setLinearVelocity (m_cube, Vector3r(sin(t), cos(t), sin(2*t)*0.5));
        sys.setLinearVelocity (m_cube, Vector3r(0,0,0));
        sys.applyForce        (m_cube, 1e10 * Vector3r(0,0,9.81), Vector3r::Zero());  // counteract gravity
    }

    private:
       entt::entity m_cube;
};
