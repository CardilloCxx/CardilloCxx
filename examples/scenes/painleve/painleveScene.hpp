#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Painlevé rod scene (moved into its own folder)
class PainleveScene : public SceneBase {
public:
    PainleveScene() = default;
    ~PainleveScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        PhysicsSystem::Cube ground;
        const real_t groundHalfThickness = (real_t)0.01;
        const real_t groundHalfSize = (real_t)50.0;
        ground.center = Vector3r(0.0, 0.0, -groundHalfThickness);
        ground.halfExtents = Vector3r(groundHalfSize, groundHalfSize, groundHalfThickness);
        ground.q = Quaternion4r::Identity();
        sys.addObstacleBody(ground);

        const real_t phi = (real_t)(31.0 * M_PI / 180.0);
        const real_t mass = (real_t)1.0;
        const real_t l = (real_t)1.0;
        const real_t thickness = (real_t)0.2;
        const Vector3r position((real_t)0.0, (real_t)0.0, l * std::sin(phi) + 1.0 * thickness * std::cos(phi));
        PhysicsSystem::Capsule rodShape;
        rodShape.radius = thickness * (real_t)0.5;
        rodShape.halfLength = l;
        const Quaternion4r capsuleOrientation(Eigen::AngleAxis<real_t>(M_PI_2 -phi, Vector3r::UnitY()));
        const real_t v0 = (real_t)30.0;
        const Vector3r linearVelocity(-v0, 0.0, 0.0);
        sys.addRigidBodyCapsule(mass, position, capsuleOrientation, linearVelocity, Vector3r::Zero(), rodShape);

        const auto& reg = sys.ecs();
        entt::entity rodEntity = entt::null;
        for (auto entity : reg.view<cardillo::PhysicsSystem::C_RigidBodyTag>()) {
            rodEntity = entity; break;
        }
        
        // Print Inertia and other stats from physics system:
        const real_t* inertiaDiag = reg.get<cardillo::PhysicsSystem::C_InertiaDiag>(rodEntity).I.data();
        PetscPrintf(PETSC_COMM_WORLD, "Painleve rod inertia diag: Ixx=%.6f, Iyy=%.6f, Izz=%.6f\n",
                    inertiaDiag[0], inertiaDiag[1], inertiaDiag[2]);
        PetscPrintf(PETSC_COMM_WORLD, " 1/3 m L^2 = %.6f\n", (real_t)(1.0/3.0) * mass * l * l);

    }
};
