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

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        // Mesh path (relative to repository root). Adjust if your runtime cwd differs.
        const std::string meshPath = "res/meshes/chain_link.obj";
        const std::string spoolMeshPath = "res/meshes/spool.obj";

        // Spool
        auto spool_shape = physics::MeshShape{spoolMeshPath, Vector3r(0.075, 0.075, 0.075)};
        auto spool_state = physics::RigidState{Vector3r(0.0, 0.0, 1.0)};
        auto spool_props = physics::RigidProps{1.0};
        m_spool = engine.addRigidBody(spool_shape, spool_state, spool_props);
        // Kinematically driven spool: fixed position, rotation with 1.5 rad/s about the drum axis
        // (world y, the y-axis of the spool mesh). The pose is taken from the created body, since a
        // dynamic mesh body is shifted to its center of mass and rotated to its principal axes (it
        // differs from spool_state); the rotation is prescribed relative to it in world axes.
        const Vector3r spoolPos = engine.ecs().get<C_Position3>(m_spool).value;
        const Quaternion4r spoolOri = engine.ecs().get<C_Orientation>(m_spool).value;
        engine.makeStatic(m_spool);
        engine.addTrajectory(m_spool, [spoolPos, spoolOri](real_t t) -> TrajectoryPose {
            return TrajectoryPose{spoolPos, Quaternion4r(Eigen::AngleAxis<real_t>((real_t)1.5 * t, Vector3r::UnitY())) * spoolOri};
        });

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

            physics::MeshShape shape{meshPath, scale};
            physics::RigidState state; state.position = pos; state.setOrientation(ori); state.linearVelocity = vlin; state.angularVelocity = omega;
            physics::RigidProps props; props.mass = mass;
            auto link = engine.addRigidBody(shape, state, props);

            if (i == 0) {
                // Hang the first link at the bottom of the spool core. The attachment points are given
                // in world coordinates and converted to the body frames, which are the principal-axes
                // frames of the meshes (they differ from the nominal mesh frames).
                const Vector3r spoolAttach = spoolPos - Vector3r(0.0, 0.0, spool_radius);  // bottom of the core
                const Vector3r linkPos = engine.ecs().get<C_Position3>(link).value;
                const Vector3r linkAttach = linkPos + Vector3r(0.0, 0.0, linkHeight / 2);  // top end of the link
                const Matrix33r R_spool = spoolOri.toRotationMatrix();
                const Matrix33r R_link = engine.ecs().get<C_Orientation>(link).value.toRotationMatrix();
                engine.addLinearDistanceConstraint(m_spool, link, R_spool.transpose() * (spoolAttach - spoolPos), R_link.transpose() * (linkAttach - linkPos));
            }
        }
    }


private:
    entt::entity m_spool{entt::null};
};
