#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Chain scene: loads a chain link mesh and creates a hanging chain of links.
// The top link is static and placed at height 1.0 m above ground. Subsequent
// links hang below it with alternating 90 degree rotations (about the X axis)
// to interlock. The mesh is assumed to represent a single chain-link of unit
// canonical height; we scale so that the link height in world space is 1.0 m.
class ChainScene : public SceneBase {
public:
    ChainScene() = default;
    ~ChainScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane (thin cube) for visual + collision
        PhysicsSystem::Cube ground;
        const real_t groundHalfThickness = (real_t)0.01;
        const real_t groundHalfSize = (real_t)50.0;
        ground.center = Vector3r(0.0, 0.0, -groundHalfThickness);
        ground.halfExtents = Vector3r(groundHalfSize, groundHalfSize, groundHalfThickness);
        ground.q = Quaternion4r::Identity();
        sys.addObstacleBody(ground);

        // Mesh path (relative to repository root). Adjust if your runtime cwd differs.
        const std::string meshPath = "res/meshes/chain_link.obj";

        // Desired physical height of one link (meters)
        const int N = 80;
        const real_t linkHeight = (real_t)0.05;
        const real_t spacing = (real_t)0.75 * linkHeight;
        const real_t startHeight = linkHeight + (N * spacing);
        const Vector3r scale(1.0* linkHeight, 1.0* linkHeight, 1.0 * linkHeight);

        // Top static link at z = 1.0 m
        const Vector3r topPos(0.0, 0.0, startHeight);
        const Quaternion4r topOri = Eigen::AngleAxis<real_t>(M_PI_2, Vector3r::UnitX()) * Quaternion4r::Identity();
        sys.addObstacleMesh(topPos, topOri, meshPath, scale);

        // Hanging dynamic links: place N links below the top link
        const real_t mass = (real_t)0.5;
        for (int i = 0; i < N; ++i) {
            Vector3r pos(0.0, spacing * (i + 1), topPos.z()); // Rotate position 90 degrees around Y
            // Alternate 0 and 90 degree rotation about X to interlock
            real_t angle = (i % 2 == 1) ? 0.0 : (real_t)(M_PI_2);
            Quaternion4r ori(Eigen::AngleAxis<real_t>(M_PI_2, Vector3r::UnitX())  * Eigen::AngleAxis<real_t>(M_PI_2, Vector3r::UnitY()) * Eigen::AngleAxis<real_t>(angle, Vector3r::UnitX()));
            // Start at rest
            Vector3r vlin = Vector3r::Zero();
            Vector3r omega = Vector3r::Zero();

            if (i == N-1) {
                // Give the last link an initial push
                vlin = Vector3r(100.0, 0.0, 100.0);
            }

            sys.addRigidBodyMesh(mass, pos, ori, vlin, omega, meshPath, scale);
        }
    }
};
