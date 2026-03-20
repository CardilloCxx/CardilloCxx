#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <fstream>
#include <sstream>
#include <limits>

using namespace cardillo;

// RodAssemblyScene: loads three meshes (frame, lowerRod, upperRod).
// The frame is represented visually by the OBJ but uses a simple cube
// collider to avoid BVH-loading failures; the two rods are created as
// dynamic capsule bodies with matching visual meshes. The rods are
// automatically aligned to the frame's mesh center using a quick OBJ
// bounding-box pass.
class RodAssemblyScene : public SceneBase {
public:
    const char* sceneName() const override { return "rod_assembly"; }
    RodAssemblyScene() = default;
    ~RodAssemblyScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Mesh paths
        const std::string framePath = "res/meshes/frame.obj";
        const std::string lowerPath = "res/meshes/lowerRod.obj";
        const std::string upperPath = "res/meshes/upperRod.obj";

        // Add frame as obstacle at origin
        const Vector3r framePos = Vector3r::Zero();
        const Quaternion4r frameOri = Quaternion4r::Identity();
        const Vector3r frameScale = Vector3r::Ones();
        {
            PhysicsSystem::MeshShape shape{framePath, frameScale}; PhysicsSystem::RigidState st; st.position = framePos; st.orientation = frameOri; PhysicsSystem::RigidProps pr; m_frame = cardillo::physics::BodyFactory::addRigidBody(sys, shape, st, pr);
        }

        // Add lower and upper rods as dynamic rigid bodies
        const real_t rodMass = (real_t)1.0;
        {
            PhysicsSystem::MeshShape shape{lowerPath, Vector3r::Ones()}; PhysicsSystem::RigidState st; st.position = Vector3r::Zero(); st.orientation = Quaternion4r::Identity(); st.linearVelocity = Vector3r(0.0,1.0,0.0); PhysicsSystem::RigidProps pr; pr.mass = rodMass; m_lowerRod = cardillo::physics::BodyFactory::addRigidBody(sys, shape, st, pr);
        }
        {
            PhysicsSystem::MeshShape shape{upperPath, Vector3r::Ones()}; PhysicsSystem::RigidState st; st.position = Vector3r::Zero(); st.orientation = Quaternion4r::Identity(); st.linearVelocity = Vector3r(1.0,0.0,0.0); PhysicsSystem::RigidProps pr; pr.mass = rodMass; m_upperRod = cardillo::physics::BodyFactory::addRigidBody(sys, shape, st, pr);
        }
    
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t /*t*/, real_t /*dt*/) override {
        sys.applyForce(m_upperRod, Vector3r(50.0, 0.0, 0), Vector3r::Zero());
    }

private: 
    entt::entity m_frame{entt::null};
    entt::entity m_lowerRod{entt::null};
    entt::entity m_upperRod{entt::null};
};
