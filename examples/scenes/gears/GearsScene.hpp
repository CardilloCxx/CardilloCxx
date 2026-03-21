#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Gears scene.
class GearsScene : public SceneBase {
public:
    const char* sceneName() const override { return "gears"; }
    GearsScene() = default;
    ~GearsScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        auto& sys = engine.world();
        using namespace cardillo;
        using namespace cardillo::misc;

        //floor 
        auto floor = engine.addStaticBody(physics::CubeShape{Vector3r{5,5,1}}, physics::RigidState{Vector3r{0,0,-1.0}});

        // gear system 1
        {
            real_t density = 5000.0;
            m_gear1 = engine.addRigidBody(physics::MeshShape("res/meshes/gears2/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, m_gear1, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));

            auto gear2 = engine.addRigidBody(physics::MeshShape("res/meshes/gears2/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, gear2, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));
        }

        // gear system 2
        {
            real_t density = 5000.0;
            auto gear3 = engine.addRigidBody(physics::MeshShape("res/meshes/gears3/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,0,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, gear3, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,1,0));

            m_gear4 = engine.addRigidBody(physics::MeshShape("res/meshes/gears3/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,0,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, m_gear4, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));
        }

        // gear system 3
        {
            real_t density = 5000.0;
            m_gear5 = engine.addRigidBody(physics::MeshShape("res/meshes/gears4/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{-1,0,0}}, physics::RigidProps::withDensity(density));
            auto rot5 = Quaternion4r(sys.getPosition(m_gear5).tail<4>()).toRotationMatrix().transpose();
            sys.addHingeConstraint(floor, m_gear5, Vector3r(0,0,0),  rot5 * Vector3r(-0.055,-0.055, 0), Vector3r(0,0,1));

            auto gear6 = engine.addRigidBody(physics::MeshShape("res/meshes/gears4/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{-1,0,0}}, physics::RigidProps::withDensity(density));
            auto rot6 = Quaternion4r(sys.getPosition(gear6).tail<4>()).toRotationMatrix();
            sys.addHingeConstraint(floor, gear6, Vector3r(0,0,0), rot6 * Vector3r(-0.05,0.05, 0), Vector3r(0,0,1));
        }

        // gear system 4
        {
            real_t density = 5000.0;
            auto gear7 = engine.addRigidBody(physics::MeshShape("res/meshes/gears5/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{0,1,0}}, physics::RigidProps::withDensity(density));
            auto rot7 = Quaternion4r(sys.getPosition(gear7).tail<4>()).toRotationMatrix();
            sys.addHingeConstraint(floor, gear7, Vector3r(0,0,0), rot7 * Vector3r(-0.015,0.07,0), Vector3r(0,0,1), 50, 1);

            m_gear8 = engine.addRigidBody(physics::MeshShape("res/meshes/gears5/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{0,1,0}}, physics::RigidProps::withDensity(density));
            auto rot8 = Quaternion4r(sys.getPosition(m_gear8).tail<4>()).toRotationMatrix();
            sys.addHingeConstraint(floor, m_gear8, Vector3r(0,0,0), rot8 * Vector3r(0,0,0), Vector3r(0,0,1));
        }

        // gear system 5
        {
            real_t density = 5000.0;
            auto gear10 = engine.addRigidBody(physics::MeshShape("res/meshes/gears6/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{0,-1,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, gear10, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));

            m_gear9 = engine.addRigidBody(physics::MeshShape("res/meshes/gears6/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{0,-1,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, m_gear9, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));
        }

        // gear system 6
        {
            real_t density = 5000.0;
            m_gear10 = engine.addRigidBody(physics::MeshShape("res/meshes/gears7/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,1,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, m_gear10, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));

            auto gear11 = engine.addRigidBody(physics::MeshShape("res/meshes/gears7/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,1,0}}, physics::RigidProps::withDensity(density));
            sys.addHingeConstraint(floor, gear11, Vector3r(0,0,0), Vector3r(0,0,0), Vector3r(0,0,1));
        }
        
        // gear system 7
        {
            real_t density = 5000.0;
            m_gear12 = engine.addRigidBody(physics::MeshShape("res/meshes/gears8/gear1.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,-1,0}}, physics::RigidProps::withDensity(density));
            auto rot12 = Quaternion4r(sys.getPosition(m_gear12).tail<4>()).toRotationMatrix().transpose();
            sys.addHingeConstraint(floor, m_gear12, Vector3r(0,0,0), rot12 * Vector3r(0.022,-0.012,0), Vector3r(0,0,1));

            auto gear13 = engine.addRigidBody(physics::MeshShape("res/meshes/gears8/gear2.obj", Vector3r(0.01, 0.01, 0.01)),physics::RigidState{Vector3r{1,-1,0}}, physics::RigidProps::withDensity(density));
            auto rot13 = Quaternion4r(sys.getPosition(gear13).tail<4>()).toRotationMatrix().transpose();
            sys.addHingeConstraint(floor, gear13, Vector3r(0,0,0), rot13 * Vector3r(-0.019,0.01,0), Vector3r(0,0,1));
        }
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        auto& sys = engine.world();
        sys.applyForce(m_gear1, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear4, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear5, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear8, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear9, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear10, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 5.0));
        sys.applyForce(m_gear12, Vector3r(0.0, 0.0, 0.0), Vector3r(0.0, 0.0, 10.0));
    }

    private:
        entt::entity m_gear1 = {entt::null};
        entt::entity m_gear4 = {entt::null};
        entt::entity m_gear5 = {entt::null};
        entt::entity m_gear8 = {entt::null};
        entt::entity m_gear9 = {entt::null};
        entt::entity m_gear10 = {entt::null};
        entt::entity m_gear12 = {entt::null};
};
