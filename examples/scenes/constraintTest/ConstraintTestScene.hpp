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
        auto floor = sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(15.0, 15.0, 0.1)), PhysicsSystem::RigidState(Vector3r(0,0,-0.1)));
        const real_t cubeHalfExt = 0.5;

        // Cube with fixed constraint
        m_cube = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(cubeHalfExt, cubeHalfExt, cubeHalfExt)),
            PhysicsSystem::RigidState(Vector3r(0,0,4), Vector3r(0,0,0), Quaternion4r::UnitRandom(),
                 Vector3r(1, 0.2, 0.3)),
            // PhysicsSystem::RigidProps::withDensity(2500.0)
            PhysicsSystem::RigidProps(1e10) 
        );

        // small cube allowed to move in a plane
        auto cube2 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, cubeHalfExt + 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, cube2);
        sys.disableCollisionBetween(m_cube, cube2);

        // small cube allowed to move in a plane
        auto cube3 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(-cubeHalfExt - 0.05, cubeHalfExt + 0.05, cubeHalfExt + 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        Vector3r K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(2) = 1000;
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube3, physics::JointFrame(cube3), K_trans);
        sys.disableCollisionBetween(m_cube, cube3);

        auto cube4 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(-cubeHalfExt - 0.05, -cubeHalfExt - 0.05, cubeHalfExt + 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(1) = 1000;
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube4, physics::JointFrame(cube4), K_trans);
        sys.disableCollisionBetween(m_cube, cube4);

        auto cube5 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, -cubeHalfExt - 0.05, cubeHalfExt + 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(0) = 1000;
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube5, physics::JointFrame(cube5), K_trans);
        sys.disableCollisionBetween(m_cube, cube5);

        auto cube6 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, -cubeHalfExt - 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube6, 
            physics::JointFrame(Vector3r(cubeHalfExt, cubeHalfExt, -cubeHalfExt), Matrix33r::Identity(), m_cube), K_trans);


        auto plate = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.3,0.3,0.02)),
            PhysicsSystem::RigidState(Vector3r(0,0,cubeHalfExt + 0.02 + 0.05), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::HingeConstraint>(sys.ecs(), m_cube, plate, 
            physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), plate)
        );

        auto cube7 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(0.2,0.2,0.05), plate, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), plate, cube7);
        sys.disableCollisionBetween(plate, cube7);
        
        
        const real_t wallHeight = 0.05;
        const real_t wallLength = 0.45;
        const real_t wallThickness = 0.01;

        auto wall1 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(wallHeight, wallLength, wallThickness)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + wallHeight, 0, wallLength), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, wall1);
        sys.disableCollisionBetween(m_cube, wall1); 

        auto wall2 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(wallHeight, wallLength, wallThickness)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + wallHeight, 0, -wallLength), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, wall2);
        sys.disableCollisionBetween(m_cube, wall2);

        auto wall3 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(wallHeight, wallThickness, wallLength)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + wallHeight, wallLength, 0), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, wall3);
        sys.disableCollisionBetween(m_cube, wall3);

        auto wall4 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(wallHeight, wallThickness, wallLength)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05, -wallLength, 0), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_cube, wall4);
        sys.disableCollisionBetween(m_cube, wall4);

        auto cube8 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.05,0.05,0.05)),
            PhysicsSystem::RigidState(Vector3r(cubeHalfExt + 0.05 + 0.01, 0, 0), m_cube, sys.ecs()),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_cube, cube8, physics::JointFrame(cube8),
            Vector3r(std::numeric_limits<real_t>::infinity(), 100, 100)
        );
        sys.addConstraint<physics::RotationConstraint>(sys.ecs(), m_cube, cube8, physics::JointFrame(cube8));
        sys.disableCollisionBetween(m_cube, cube8);

        // -----------------------------------------------------------------
        // Double pendulum built from three spheres
        // -----------------------------------------------------------------

        const real_t sphereRadius = 0.00005;

        // Top sphere: static, acts as fixed anchor in world
        auto topSphere = sys.addStaticBody(
            PhysicsSystem::SphereShape(sphereRadius),
            PhysicsSystem::RigidState(Vector3r(-1.0, 0.0, 2.5))
        );
        sys.makeStatic(topSphere);

        // First pendulum mass
        auto midSphere = sys.addRigidBody(
            PhysicsSystem::SphereShape(sphereRadius),
            PhysicsSystem::RigidState(Vector3r(-1.0, 0.0, 2.2), Vector3r(1.0, 0.0, 0.0)),
            PhysicsSystem::RigidProps(1.0)
        );

        // Second pendulum mass
        auto bottomSphere = sys.addRigidBody(
            PhysicsSystem::SphereShape(sphereRadius),
            PhysicsSystem::RigidState(Vector3r(-1.0, 0.0, 1.9), Vector3r(-4.0, 0.0, 0.0)),
            PhysicsSystem::RigidProps(1.0)
        );
        sys.track(midSphere, "mid_sphere");
        sys.track(bottomSphere, "bottom_sphere");

        // Hinge between top (static) sphere and first mass
        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(), topSphere, midSphere,
            physics::JointFrame::fromAxis(
                Vector3r(0,0,0), // hinge between top and mid
                Vector3r::UnitY(),          // rotation axis
                topSphere
            )
        );

        // Hinge between first and second mass
        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(), midSphere, bottomSphere,
            physics::JointFrame::fromAxis(
                Vector3r(0,0,0), // hinge between mid and bottom
                Vector3r::UnitY(),
                midSphere
            )
        );

        sys.disableCollisionBetween(topSphere, midSphere);
        sys.disableCollisionBetween(midSphere, bottomSphere);
        sys.disableCollisionBetween(topSphere, bottomSphere);

        // -------------------------------------------------------------
        // Two levers hinged to the ground with different stiffnesses
        // -------------------------------------------------------------

        // First lever: stiffer hinge
        m_lever_1 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            PhysicsSystem::RigidState(Vector3r(-1.5, 1.0, 0.2)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        // Hinge at the left end of the lever, anchored in the world frame
        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(),floor,  m_lever_1,
            physics::JointFrame::fromAxis(
                Vector3r(-0.4, 0.0, 0.0),   // world anchor position
                Vector3r::UnitZ(),          // hinge axis
                m_lever_1
            ),
            1e0
        );

        // Second lever: softer hinge
        m_lever_2 = sys.addRigidBody(
            PhysicsSystem::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            PhysicsSystem::RigidState(Vector3r(-1.5, -1.0, 0.2)),
            PhysicsSystem::RigidProps::withDensity(2500.0)
        );

        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(), floor, m_lever_2,
            physics::JointFrame::fromAxis(
                Vector3r(-0.4, 0.0, 0.0),
                Vector3r::UnitZ(),
                m_lever_2
            ),
            2e0
        );
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        
        const Vector3r torque(0.0, 0.0, 5.0);

        sys.applyForce(m_lever_1, Vector3r::Zero(), torque);
        sys.applyForce(m_lever_2, Vector3r::Zero(), torque);

        sys.applyForce(m_cube, 1e10 * Vector3r(0,0,9.81), Vector3r::Zero());  // counteract gravity
    }

    private:
       entt::entity m_cube;
       entt::entity m_lever_1;
       entt::entity m_lever_2;
};
