#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Gears scene.
class PendulumScene : public SceneBase {
public:
    const char* sceneName() const override { return "pendulum"; }
    PendulumScene() = default;
    ~PendulumScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // floor 
        auto floor = sys.addStaticBody(
            PhysicsSystem::CubeShape{Vector3r{5,5,0.1}}, 
            PhysicsSystem::RigidState{Vector3r{0,0,-3.5}}
        );

        // capsule 1
        const real_t phi1 = M_PI / 2;
        const real_t mass1 = 1.0;
        const real_t l1 = 1.0;
        const real_t thickness1 = 0.2;
        PhysicsSystem::CapsuleShape shape1{thickness1 * 0.5, l1};
        
        PhysicsSystem::RigidState state1; 
        state1.position = Vector3r(l1, 0.0, 0.0);
        state1.orientation = Quaternion4r(AngleAxis3r(phi1, Vector3r::UnitY()));
        state1.linearVelocity = Vector3r::Zero(); 
        state1.angularVelocity = Vector3r::Zero();

        PhysicsSystem::RigidProps props1;
        props1.mass = mass1;

        auto capsule1 = sys.addRigidBody(shape1, state1, props1);

        // capsule 2
        const real_t phi2 = M_PI / 2;
        const real_t mass2 = 1.0;
        const real_t l2 = 1.0;
        const real_t thickness2 = 0.1;
        PhysicsSystem::CapsuleShape shape2{thickness2 * 0.5, l2};
        
        PhysicsSystem::RigidState state2; 
        state2.position = Vector3r(2 * l1 + l2, 0.0, 0.0);
        state2.orientation = Quaternion4r(AngleAxis3r(phi2, Vector3r::UnitY()));
        state2.linearVelocity = Vector3r::Zero(); 
        state2.angularVelocity = Vector3r::Zero();

        PhysicsSystem::RigidProps props2;
        props2.mass = mass2;

        auto capsule2 = sys.addRigidBody(shape2, state2, props2);

        // disable collisions
        // sys.disableCollisionBetween(floor, capsule1);
        // sys.disableCollisionBetween(floor, capsule2);
        sys.disableCollisionBetween(capsule1, capsule2);

        // hinge origin <> capsule1
        // sys.addConstraint<physics::RigidConstraint>(sys.ecs(), floor, capsule1, Vector3r::Zero());
        Vector3r k_axis1 = Vector3r::Constant(std::numeric_limits<real_t>::max());
        k_axis1(2) = 1e-1; // small stiffness on y-axis
        Vector3r d_axis1 = Vector3r::Zero(); // no damping
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), floor, capsule1, Vector3r::Zero(), std::nullopt, k_axis1, d_axis1);

        // hinge capsule1 <> capsule2
        // sys.addConstraint<physics::RigidConstraint>(sys.ecs(), capsule1, capsule2, Vector3r(2 * l1,0,0));
        Vector3r k_axis2 = Vector3r::Constant(std::numeric_limits<real_t>::max());
        k_axis2(2) = 1e2; // moderate stiffness on y-axis
        Vector3r d_axis2 = Vector3r::Zero();
        // d_axis2(2) = 1e-3; // small stiffness on y-axis
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), capsule1, capsule2, Vector3r(2 * l1,0,0), std::nullopt, k_axis2, d_axis2);
    }

    // void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
    //     sys.applyForce(m_gear1, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear4, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear5, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear8, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear9, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear10, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
    //     sys.applyForce(m_gear12, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 10.0));
    // }

    // private:
    //     entt::entity m_gear1 = {entt::null};
    //     entt::entity m_gear4 = {entt::null};
    //     entt::entity m_gear5 = {entt::null};
    //     entt::entity m_gear8 = {entt::null};
    //     entt::entity m_gear9 = {entt::null};
    //     entt::entity m_gear10 = {entt::null};
    //     entt::entity m_gear12 = {entt::null};
};
