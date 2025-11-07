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

    // Ground plane (static)
    PhysicsSystem::PlaneShape groundShape{Vector3r(0,0,1), Vector3r(0,1,0), (real_t)6.0, (real_t)6.0};
    PhysicsSystem::RigidState groundState; groundState.position = Vector3r(0,0,0); groundState.orientation = Quaternion4r::Identity();
    sys.addStaticBody(groundShape, groundState);

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
            PhysicsSystem::CubeShape wall1Shape{Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0)};
            PhysicsSystem::RigidState wall1State; wall1State.position = Vector3r((rows * dia) / 2.0 - r, -wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
            wall1State.orientation = Quaternion4r(Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ()));
            PhysicsSystem::RigidProps wall1Props; sys.addRigidBody(wall1Shape, wall1State, wall1Props);

            PhysicsSystem::CubeShape wall2Shape{Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0)};
            Vector3r wall2Center = Vector3r((rows * dia) / 4.0 - r, (rows * dia * std::sqrt(3.0)) / 4.0 + wall_thickness / 2.0 + r, (dia * std::sqrt(6.0)) / 6.0);
            wall2Center -= Vector3r(std::cos(M_PI / 3.0) * 1.25 * r, std::sin(M_PI / 3.0) * 1.25 * r, 0.0);
            PhysicsSystem::RigidState wall2State; wall2State.position = wall2Center; wall2State.orientation = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI / 3.0, Vector3r::UnitZ()));
            PhysicsSystem::RigidProps wall2Props; sys.addRigidBody(wall2Shape, wall2State, wall2Props);

            PhysicsSystem::CubeShape wall3Shape{Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0)};
            Vector3r wall3Center = Vector3r((rows * dia) * 3.0 / 4.0 - r / 2.0, (rows * dia * std::sqrt(3.0)) / 4.0 - wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
            wall3Center += Vector3r(std::cos(M_PI / 3.0) * r * 1.25, std::sin(M_PI / 3.0) * r * 1.25, 0.0);
            PhysicsSystem::RigidState wall3State; wall3State.position = wall3Center; wall3State.orientation = Quaternion4r(Eigen::AngleAxis<real_t>(-M_PI / 3.0, Vector3r::UnitZ()));
            PhysicsSystem::RigidProps wall3Props; sys.addRigidBody(wall3Shape, wall3State, wall3Props);
        }

        // Add a heavy point mass above to drop into pyramid
        sys.addPointMass(10000.0, base_center + Vector3r(r * rows, rows * r / std::sqrt(3.0), 5.0), Vector3r(0.0, 0.0, -10.0), 10*r);
    }
};
