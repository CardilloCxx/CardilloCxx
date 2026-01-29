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

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        real_t density = 7750.0; // kg/m^3
        Vector3r scale{1.0, 1.0, 1.0};

        // base link
        m_base = sys.addStaticBody(
            PhysicsSystem::MeshShape("res/meshes/double_pendulum/base_link.stl", scale), 
            PhysicsSystem::RigidState{}
        );

        // first rigid body
        m_rb1 = sys.addRigidBody(
            PhysicsSystem::MeshShape("res/meshes/double_pendulum/link1.stl", scale),
            // PhysicsSystem::RigidState{Vector3r{3.0299e-03, 2.6279e-13, 2.9120e-02}}, // TODO: What is the initial configuration?
            PhysicsSystem::RigidState{
                Vector3r{0.0, 0.0, 0.035}, // initial position
                Quaternion4r{Eigen::AngleAxis<real_t>(-M_PI / 2, Vector3r::UnitX())}, // initial orientation
            },
            PhysicsSystem::RigidProps::withDensity(density)
        );

        // no collision between base and first rigid body
        sys.disableCollisionBetween(m_base, m_rb1);

        // add revolute joint between base and first rigid body
        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(), 
            m_base, 
            m_rb1, 
            physics::JointFrame(m_base)
        );

        // second rigid body
        m_rb2 = sys.addRigidBody(
            PhysicsSystem::MeshShape("res/meshes/double_pendulum/link2.stl", scale),
            PhysicsSystem::RigidState{
                Vector3r{0.0, 0.1, 0.04}, // initial position
                Quaternion4r{Eigen::AngleAxis<real_t>(-M_PI / 2, Vector3r::UnitX())}, // initial orientation
            },
            PhysicsSystem::RigidProps::withDensity(density)
        );

        // no collision between first and second rigid body
        sys.disableCollisionBetween(m_rb1, m_rb2);

        // add revolute joint between base and first rigid body
        sys.addConstraint<physics::HingeConstraint>(
            sys.ecs(), 
            m_rb1, 
            m_rb2, 
            physics::JointFrame(m_rb2)
        );
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        // nothing done here
    }

    private:
        entt::entity m_base = {entt::null};
        entt::entity m_rb1 = {entt::null};
        entt::entity m_rb2 = {entt::null};
};
