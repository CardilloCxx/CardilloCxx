#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Chain scene: loads a chain link mesh and creates a hanging chain of links.
// The top link is static and placed at height 1.0 m above ground.
class ChainScene : public SceneBase {
public:
    const char* sceneName() const override { return "chain"; }
    ChainScene() = default;
    ~ChainScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Mesh path (relative to repository root). Adjust if your runtime cwd differs.
        const std::string meshPath = "res/meshes/chain_link.obj";
        const std::string spoolMeshPath = "res/meshes/spool.obj";

        // Spool
        auto spool_shape = PhysicsSystem::MeshShape{spoolMeshPath, Vector3r(0.075, 0.075, 0.075)};
        auto spool_state = PhysicsSystem::RigidState{Vector3r(0.0, 0.0, 1.0)};
        auto spool_props = PhysicsSystem::RigidProps{1e10}; // heavy to stay put
        m_spool = sys.addRigidBody(spool_shape, spool_state, spool_props);

        // Desired physical height of one link (meters)
        const int N = 150;
        const real_t linkHeight = (real_t)0.05;
        const real_t spacing = (real_t)0.75 * linkHeight;
        const real_t spool_radius = 0.075;
        const real_t startHeight = 1.0 - spool_radius;
        const Vector3r scale(1.0* linkHeight, 1.0* linkHeight, 1.0 * linkHeight);

        // Hanging dynamic links: place N links below the top link

        const real_t mass = (real_t)0.5;
        for (int i = 0; i < N; ++i) {
            Vector3r pos(0.0, 0.0, startHeight - spacing * (i + 1)); // Rotate position 90 degrees around Y
            // Alternate 0 and 90 degree rotation about X to interlock
            real_t angle = (i % 2 == 0) ? 0.0 : (real_t)(M_PI_2);
            Quaternion4r ori(Eigen::AngleAxis<real_t>(angle, Vector3r::UnitZ()));
            // Start at rest
            Vector3r vlin = Vector3r::Zero();
            Vector3r omega = Vector3r::Zero();

            PhysicsSystem::MeshShape shape{meshPath, scale};
            PhysicsSystem::RigidState state; state.position = pos; state.orientation = ori; state.linearVelocity = vlin; state.angularVelocity = omega;
            PhysicsSystem::RigidProps props; props.mass = mass;
            auto link = sys.addRigidBody(shape, state, props);

            if (i == 0) {
                  sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), m_spool, link, Vector3r(0.0, 0.0, -spool_radius), Vector3r(linkHeight / 2,  0.0, 0.0));
            }
        }
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        sys.setAngularVelocity(m_spool, Vector3r(1.5, 0.0, 0.0)); 
        sys.setLinearVelocity (m_spool, Vector3r::Zero());
        sys.applyForce        (m_spool, 1e10 * Vector3r(0,0,9.81), Vector3r::Zero());  // counteract gravity
    }

private:
    entt::entity m_spool{entt::null};
};
