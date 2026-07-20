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

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        // Floor (static)
        auto floor = engine.addStaticBody(physics::CubeShape(Vector3r(15.0, 15.0, 0.1)), physics::RigidState(Vector3r(0,0,-0.1)));
        const real_t cubeHalfExt = 0.5;

        // Cube with fixed constraint
        m_cube = engine.addRigidBody(            physics::CubeShape(Vector3r(cubeHalfExt, cubeHalfExt, cubeHalfExt)),
            physics::RigidState(Vector3r(0,0,4), Vector3r(0,0,0), Quaternion4r::UnitRandom(),
                 Vector3r(1, 0.2, 0.3)),
            // physics::RigidProps::withDensity(2500.0)
            physics::RigidProps(1e10) 
        );

        // small cube allowed to move in a plane
        auto cube2 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, cubeHalfExt + 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );

        engine.addRigidConstraint(m_cube, cube2);
        engine.disableCollisionBetween(m_cube, cube2);

        // small cube allowed to move in a plane
        auto cube3 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(-cubeHalfExt - 0.05, cubeHalfExt + 0.05, cubeHalfExt + 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );

        Vector3r K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(2) = 1000;
        engine.addTranslationalConstraint(m_cube, cube3, physics::JointFrame(cube3), K_trans);
        engine.disableCollisionBetween(m_cube, cube3);

        auto cube4 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(-cubeHalfExt - 0.05, -cubeHalfExt - 0.05, cubeHalfExt + 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(1) = 1000;
        engine.addTranslationalConstraint(m_cube, cube4, physics::JointFrame(cube4), K_trans);
        engine.disableCollisionBetween(m_cube, cube4);

        auto cube5 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(cubeHalfExt + 0.05, -cubeHalfExt - 0.05, cubeHalfExt + 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        K_trans(0) = 1000;
        engine.addTranslationalConstraint(m_cube, cube5, physics::JointFrame(cube5), K_trans);
        engine.disableCollisionBetween(m_cube, cube5);

        auto cube6 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(cubeHalfExt + 0.05, cubeHalfExt + 0.05, -cubeHalfExt - 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );

        K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
        engine.addTranslationalConstraint(m_cube, cube6, 
            physics::JointFrame(Vector3r(cubeHalfExt, cubeHalfExt, -cubeHalfExt), Matrix33r::Identity(), m_cube), K_trans);


        auto plate = engine.addRigidBody(            physics::CubeShape(Vector3r(0.3,0.3,0.02)),
            physics::RigidState(Vector3r(0,0,cubeHalfExt + 0.02 + 0.05), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addHingeConstraint(m_cube, plate, 
            physics::JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitZ(), plate)
        );

        auto cube7 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(0.2,0.2,0.05), plate, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addRigidConstraint(plate, cube7);
        engine.disableCollisionBetween(plate, cube7);
        
        
        const real_t wallHeight = 0.05;
        const real_t wallLength = 0.45;
        const real_t wallThickness = 0.01;

        auto wall1 = engine.addRigidBody(            physics::CubeShape(Vector3r(wallHeight, wallLength, wallThickness)),
            physics::RigidState(Vector3r(cubeHalfExt + wallHeight, 0, wallLength), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addRigidConstraint(m_cube, wall1);
        engine.disableCollisionBetween(m_cube, wall1); 

        auto wall2 = engine.addRigidBody(            physics::CubeShape(Vector3r(wallHeight, wallLength, wallThickness)),
            physics::RigidState(Vector3r(cubeHalfExt + wallHeight, 0, -wallLength), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addRigidConstraint(m_cube, wall2);
        engine.disableCollisionBetween(m_cube, wall2);

        auto wall3 = engine.addRigidBody(            physics::CubeShape(Vector3r(wallHeight, wallThickness, wallLength)),
            physics::RigidState(Vector3r(cubeHalfExt + wallHeight, wallLength, 0), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addRigidConstraint(m_cube, wall3);
        engine.disableCollisionBetween(m_cube, wall3);

        auto wall4 = engine.addRigidBody(            physics::CubeShape(Vector3r(wallHeight, wallThickness, wallLength)),
            physics::RigidState(Vector3r(cubeHalfExt + 0.05, -wallLength, 0), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addRigidConstraint(m_cube, wall4);
        engine.disableCollisionBetween(m_cube, wall4);

        auto cube8 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.05,0.05,0.05)),
            physics::RigidState(Vector3r(cubeHalfExt + 0.05 + 0.01, 0, 0), m_cube, engine.ecs()),
            physics::RigidProps::withDensity(2500.0)
        );
        engine.addTranslationalConstraint(m_cube, cube8, physics::JointFrame(cube8),
            Vector3r(std::numeric_limits<real_t>::infinity(), 100, 100)
        );
        engine.addRotationConstraint(m_cube, cube8, physics::JointFrame(cube8));
        engine.disableCollisionBetween(m_cube, cube8);

        // -----------------------------------------------------------------
        // Double pendulum built from three spheres
        // -----------------------------------------------------------------

        const real_t sphereRadius = 0.00005;

        // Top sphere: static, acts as fixed anchor in world
        auto topSphere = engine.addStaticBody(            physics::SphereShape(sphereRadius),
            physics::RigidState(Vector3r(-1.0, 0.0, 2.5))
        );
        engine.makeStatic(topSphere);

        // First pendulum mass
        auto midSphere = engine.addRigidBody(            physics::SphereShape(sphereRadius),
            physics::RigidState(Vector3r(-1.0, 0.0, 2.2), Vector3r(1.0, 0.0, 0.0)),
            physics::RigidProps(1.0)
        );

        // Second pendulum mass
        auto bottomSphere = engine.addRigidBody(            physics::SphereShape(sphereRadius),
            physics::RigidState(Vector3r(-1.0, 0.0, 1.9), Vector3r(-4.0, 0.0, 0.0)),
            physics::RigidProps(1.0)
        );

        // Hinge between top (static) sphere and first mass
        engine.addHingeConstraint(topSphere, midSphere,
            physics::JointFrame::fromAxis(
                Vector3r(0,0,0), // hinge between top and mid
                Vector3r::UnitY(),          // rotation axis
                topSphere
            )
        );

        // Hinge between first and second mass
        engine.addHingeConstraint(midSphere, bottomSphere,
            physics::JointFrame::fromAxis(
                Vector3r(0,0,0), // hinge between mid and bottom
                Vector3r::UnitY(),
                midSphere
            )
        );

        engine.disableCollisionBetween(topSphere, midSphere);
        engine.disableCollisionBetween(midSphere, bottomSphere);
        engine.disableCollisionBetween(topSphere, bottomSphere);

        // -------------------------------------------------------------
        // Two levers hinged to the ground with different stiffnesses
        // -------------------------------------------------------------

        // First lever: stiffer hinge
        m_lever_1 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            physics::RigidState(Vector3r(-1.5, 1.0, 0.2)),
            physics::RigidProps::withDensity(2500.0)
        );
        real_t mass = engine.getMass(m_lever_1)(0,0);
        real_t Iz_com = engine.getInertiaDiag(m_lever_1)(2);
        real_t Iz_hinge_1 = Iz_com + mass * 0.4*0.4; // parallel axis theorem

        // Hinge at the left end of the lever, anchored in the world frame
        engine.addHingeConstraint(floor,  m_lever_1,
            physics::JointFrame::fromAxis(
                Vector3r(-0.4, 0.0, 0.0),   // world anchor position
                Vector3r::UnitZ(),          // hinge axis
                m_lever_1
            ),
            1e0
        );

        // Second lever: softer hinge
        m_lever_2 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            physics::RigidState(Vector3r(-1.5, -1.0, 0.2)),
            physics::RigidProps::withDensity(2500.0)
        );

        mass = engine.getMass(m_lever_2)(0,0);
        Iz_com = engine.getInertiaDiag(m_lever_2)(2);
        real_t Iz_hinge_2 = Iz_com + mass * 0.4*0.4; // parallel axis theorem

        engine.addHingeConstraint(floor, m_lever_2,
            physics::JointFrame::fromAxis(
                Vector3r(-0.4, 0.0, 0.0),
                Vector3r::UnitZ(),
                m_lever_2
            ),
            2e0
        );

        // engine.track(m_lever_1, "lever_Iz=" + std::to_string(Iz_hinge_1) + "_k_=1");
        // engine.track(m_lever_2, "lever_Iz=" + std::to_string(Iz_hinge_2) + "_k_=2");

                // First lever: stiffer hinge
        m_lever_3 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            physics::RigidState(Vector3r(-1.5, 3.0, 0.2)),
            physics::RigidProps::withDensity(2500.0)
        );
        mass = engine.getMass(m_lever_3)(0,0);
        Iz_com = engine.getInertiaDiag(m_lever_3)(2);
        real_t Iz_hinge_3 = Iz_com; // parallel axis theorem

        // Hinge at the left end of the lever, anchored in the world frame
        engine.addHingeConstraint(floor,  m_lever_3,
            physics::JointFrame::fromAxis(
                Vector3r(-0.0, 0.0, 0.0),   // world anchor position
                Vector3r::UnitZ(),          // hinge axis
                m_lever_3
            ),
            1e0
        );

        // Second lever: softer hinge
        m_lever_4 = engine.addRigidBody(            physics::CubeShape(Vector3r(0.4, 0.05, 0.02)),
            physics::RigidState(Vector3r(-1.5, -3.0, 0.2)),
            physics::RigidProps::withDensity(2500.0)
        );

        mass = engine.getMass(m_lever_4)(0,0);
        Iz_com = engine.getInertiaDiag(m_lever_4)(2);
        real_t Iz_hinge_4 = Iz_com; // parallel axis theorem

        engine.addHingeConstraint(floor, m_lever_4,
            physics::JointFrame::fromAxis(
                Vector3r(-0.0, 0.0, 0.0),
                Vector3r::UnitZ(),
                m_lever_4
            ),
            2e0
        );

        // engine.track(m_lever_3, "lever_Iz=" + std::to_string(Iz_hinge_3) + "_k_=1");
        // engine.track(m_lever_4, "lever_Iz=" + std::to_string(Iz_hinge_4) + "_k_=2");
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        const Vector3r torque(0.0, 0.0, 5.0);

        engine.applyForce(m_lever_1, Vector3r::Zero(), torque);
        engine.applyForce(m_lever_2, Vector3r::Zero(), torque);
        engine.applyForce(m_lever_3, Vector3r::Zero(), torque);
        engine.applyForce(m_lever_4, Vector3r::Zero(), torque);

        engine.applyForce(m_cube, 1e10 * Vector3r(0,0,9.81), Vector3r::Zero());  // counteract gravity
    }

    private:
        entt::entity m_cube;
        entt::entity m_lever_1;
        entt::entity m_lever_2;
        entt::entity m_lever_3;
        entt::entity m_lever_4;
};
