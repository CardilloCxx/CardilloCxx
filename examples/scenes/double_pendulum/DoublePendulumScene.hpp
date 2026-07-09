#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Gears scene.
class DoublePendulumScene : public SceneBase {
public:
    const char* sceneName() const override { return "double_pendulum"; }
    DoublePendulumScene() = default;
    ~DoublePendulumScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;
        using namespace cardillo::physics;

        real_t density = 7750.0; // kg/m^3
        Vector3r scale{1.0, 1.0, 1.0};

        // base link
        m_base = engine.addStaticBody(
            MeshShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/double_pendulum/base_link.stl", scale),
            RigidState{}
        );

        // first rigid body
        m_rb1 = engine.addRigidBody(
            MeshShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/double_pendulum/link1.stl", scale),
            // When using meshes the actually position and orientation differ 
            // to values passed here, as the mesh COM and principal axes are 
            // used. The passed values are passsed are offset from the model 
            // origin, when passing zeroes the model appears exactly where it 
            // is in the model, yet the COM may not be zero. 
            RigidState{
                Vector3r{0.007, 0.0, 0.035}, // initial position
                Quaternion4r{AngleAxis3r(-M_PI / 2, Vector3r::UnitX())}, // initial orientation
            },
            RigidProps::withDensity(density)
        );

        // no collision between base and first rigid body
        engine.disableCollisionBetween(m_base, m_rb1);

        // add revolute joint between base and first rigid body
        engine.addHingeConstraint(m_base, 
            m_rb1, 
            physics::JointFrame::fromAxis(
                Vector3r(0.02, 0.0002, 0.0349), // World position of the joint
                Vector3r::UnitX() // World axis of the joint
            ),
            (real_t)0.0, // optional: axis spring stiffness
            (real_t)0.0  // optional: axis damper damping
        );

        // second rigid body
        m_rb2 = engine.addRigidBody(
            MeshShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/double_pendulum/link2.stl", scale),
            RigidState{
                Vector3r{0.007, 0.1, 0.04 - 0.004}, // initial position
                Quaternion4r{Eigen::AngleAxis<real_t>(-M_PI / 2, Vector3r::UnitX())}, // initial orientation
            },
            RigidProps::withDensity(density)
        );

        // no collision between first and second rigid body
        engine.disableCollisionBetween(m_rb1, m_rb2);

        // add revolute joint between base and first rigid body
        engine.addHingeConstraint(m_rb1, 
            m_rb2, 
            physics::JointFrame::fromAxis(
                Vector3r(0.0058, 0.0991, 0.0344), // World position of the joint
                Vector3r::UnitX() // World axis of the joint
            ),
            (real_t)0.0, // optional: axis spring stiffness
            (real_t)0.0  // optional: axis damper damping
        );

        // // optionally disable collision between second rigid body and base
        // engine.disableCollisionBetween(m_base, m_rb2);
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        (void)engine;
        (void)t;
        // nothing done here
    }

    private:
        entt::entity m_base = {entt::null};
        entt::entity m_rb1 = {entt::null};
        entt::entity m_rb2 = {entt::null};
};
