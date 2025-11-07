#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

class HeightmapScene : public SceneBase {
public:
    HeightmapScene() = default;
    ~HeightmapScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        const std::string exrPath = "./res/heightmaps/mountain_height.exr";
        const real_t x_dim = (real_t)50.0;
        const real_t y_dim = (real_t)50.0;
        const real_t z_scale = (real_t)200.0;
        const real_t min_height = (real_t)0.0;
        const Vector3r pos(0.0, 0.0, 0.0);
        const Quaternion4r q = Quaternion4r::Identity();
        sys.addObstacleHeightField(pos, q, exrPath, x_dim, y_dim, z_scale, min_height);

        // Drop a few rigid-body spheres onto the mountain
        {
            const real_t radius = (real_t)0.01;
            const real_t mass = (real_t)1.0;
            const int n = 5;
            for (int i = -n; i <= n; ++i) {
                for (int j = -n; j <= n; ++j) {
                    Vector3r start((real_t)i * radius * 2.0, (real_t)j * radius * 2.0, (real_t)9.0);
                    PhysicsSystem::SphereShape shape{radius};
                    PhysicsSystem::RigidState state; state.position = start; state.orientation = Quaternion4r::Identity();
                    PhysicsSystem::RigidProps props; props.mass = mass;
                    sys.addRigidBody(shape, state, props);
                }
            }
        }
    }
};
