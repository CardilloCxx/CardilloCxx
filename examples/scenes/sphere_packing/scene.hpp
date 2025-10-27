#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// SpherePacking scene (from collisions_demo.cpp logic)
class SpherePackingScene : public SceneBase {
public:
    SpherePackingScene() = default;
    ~SpherePackingScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane
        PhysicsSystem::Plane ground;
        ground.center = Vector3r(0,0,0);
        ground.normal = Vector3r(0,0,1);
        ground.up = Vector3r(0,1,0);
        ground.sizeX = 6.0; ground.sizeY = 6.0;
        sys.addObstacleBody(ground);

        const real_t r = 0.075;
        const real_t dia = 2.0 * r;
        const int rows = 8; // reduced for quick runs
        const Vector3r base_center(0.0, 0.0, 0.0 + r);

        for (int k = 0; k < rows; ++k) {
            int n = rows - k;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n - i; ++j) {
                    real_t x = (2.0 * i + (( j + k))) * r;
                    real_t y = std::sqrt(3.0) * (j + 1.0 / 3.0 * (k)) * r;
                    real_t z = std::sqrt(6.0) * 2.0 / 3.0 * k * r;
                    Vector3r pos = base_center + Vector3r(x, y, z);
                    Vector3r vel = Vector3r::Zero();
                    sys.addPointMass(1.0, pos, vel, r * 1.01);
                }
            }
        }

        // Add three walls around the base of the pyramid
        const real_t wall_thickness = 0.1;
        {
            PhysicsSystem::Cube wall1;
            wall1.center = Vector3r((rows * dia) / 2.0 - r, -wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
            wall1.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
            wall1.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
            sys.addObstacleBody(wall1);

            PhysicsSystem::Cube wall2;
            wall2.center = Vector3r((rows * dia) / 4.0 - r, (rows * dia * std::sqrt(3.0)) / 4.0 + wall_thickness / 2.0 + r, (dia * std::sqrt(6.0)) / 6.0);
            wall2.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
            wall2.center -= Vector3r(std::cos(M_PI / 3.0) * 1.25 * r, std::sin(M_PI / 3.0) * 1.25 * r, 0.0);
            wall2.q = Eigen::AngleAxis<real_t>(M_PI / 3.0, Vector3r::UnitZ());
            sys.addObstacleBody(wall2);

            PhysicsSystem::Cube wall3;
            wall3.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
            wall3.center = Vector3r((rows * dia) * 3.0 / 4.0 - r / 2.0, (rows * dia * std::sqrt(3.0)) / 4.0 - wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
            wall3.center += Vector3r(std::cos(M_PI / 3.0) * r * 1.25, std::sin(M_PI / 3.0) * r * 1.25, 0.0);
            wall3.q = Eigen::AngleAxis<real_t>(-M_PI / 3.0, Vector3r::UnitZ());
            sys.addObstacleBody(wall3);
        }

        // Add a heavy point mass above to drop into pyramid
        sys.addPointMass(10000.0, base_center + Vector3r(r * rows, rows * r / std::sqrt(3.0), 5.0), Vector3r(0.0, 0.0, -10.0), 10*r);
    }
};
