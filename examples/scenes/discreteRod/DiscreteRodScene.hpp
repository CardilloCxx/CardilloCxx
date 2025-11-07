#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class DiscreteRodScene : public SceneBase {
public:
    DiscreteRodScene() = default;
    ~DiscreteRodScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Floor (static)
        sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(15.0, 15.0, 0.1)), PhysicsSystem::RigidState(Vector3r(0,0,-0.1)));

        const size_t segments = 50;
        const real_t length = M_PI;
        const real_t width = 0.05 * length;
        const real_t height = width / 2;
        const real_t density = 600;
        real_t E = 5e7;
        real_t nu = 0.3;

        // Circle spline with circumference=length
        const real_t radius = length / ((real_t)2 * M_PI);
        misc::CircleSpline spline(Vector3r(0,0,4), radius, Vector3r::UnitX(), Vector3r::UnitZ());
        (void) sys.createBeam(spline, segments, height, width, density, E, nu);

        // Linear Beam
        misc::LinearSpline line(Vector3r(3,-1,4), Vector3r(3,1,4));
        auto endpoints = sys.createBeam(line, segments, width, height, density, E, nu);
        sys.makeStatic(endpoints.first);
        m_rodEnd = endpoints.second;

        // Helix Beam
        misc::HelixSpline helix(Vector3r(-3,0,4), Vector3r::UnitZ(), /* radius */ 0.5, /* pitch */ 1.0, /*turns*/ 2.0);
        auto helix_endpoints = sys.createBeam(helix, segments, width, height, density, E * 10.0, nu);
        sys.makeStatic(helix_endpoints.second);

        // Beam rope between two static cubes
        auto startCube = sys.addStaticBody(PhysicsSystem::CubeShape(Vector3r(0.1,0.1,0.1)), PhysicsSystem::RigidState(Vector3r(6,-0.1,4)));
        auto endCube = sys.addRigidBody(PhysicsSystem::CubeShape(Vector3r(0.1,0.1,0.1)), PhysicsSystem::RigidState(Vector3r(6,2.1,4)), PhysicsSystem::RigidProps(1e10));
        m_endCube = endCube;

        misc::LinearSpline rope1(Vector3r(6.05, 0, 4), Vector3r(6.05, 2, 4));
        endpoints = sys.createBeam(rope1, 150, 0.03, 0.03, density, 100, 0.3, 1e7, 1e7, 0, 0, 0, Vector3r::Constant(1), Vector3r::Constant(1));
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endpoints.first, startCube);
        sys.disableCollisionBetween(endpoints.first, startCube);
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), endCube, endpoints.second);
        sys.disableCollisionBetween(endpoints.second, endCube);

        misc::LinearSpline rope2(Vector3r(5.95, 0, 4), Vector3r(5.95, 2, 4));
        endpoints = sys.createBeam(rope2, 150, 0.03, 0.03, density, 100, 0.3, 1e7, 1e7, 0, 0, 0, Vector3r::Constant(1), Vector3r::Constant(1));
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), endpoints.first, startCube);
        sys.disableCollisionBetween(endpoints.first, startCube);
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), endCube, endpoints.second);
        sys.disableCollisionBetween(endpoints.second, endCube);
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        // Apply a twisting moment at the rod end
        real_t torque_magnitude = 0.05;
        sys.applyForce(m_rodEnd, Vector3r::Zero(), Vector3r(0, -torque_magnitude, 0));
        // sys.applyForce(m_ropeEnd, Vector3r::Zero(), Vector3r(0.1, 0, 0));

        if (t < 1.0) {
            sys.setLinearVelocity(m_endCube, Vector3r(0, -1, 0));
            sys.setAngularVelocity(m_endCube, Vector3r::Zero());
        }else {
            sys.setLinearVelocity(m_endCube, Vector3r::Zero());
            sys.setAngularVelocity(m_endCube, Vector3r(0, 2, 0));
        }
        sys.applyForce(m_endCube, -sys.gravity() * sys.getMass(m_endCube).diagonal(), Vector3r::Zero());
    }

    private:
        entt::entity m_rodEnd{entt::null};
        entt::entity m_ropeEnd{entt::null};
        entt::entity m_endCube{entt::null};
};
