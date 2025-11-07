#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Dzhanibekov rod scene (moved into its own folder)
class DzhanibekovScene : public SceneBase {
public:
    DzhanibekovScene() = default;
    ~DzhanibekovScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        const real_t mass = (real_t)1.0;
        const real_t a = (real_t)3.0;
        const real_t b = (real_t)2.0;
        const real_t c = (real_t)1.0;
        const real_t eps = (real_t)1e-3;
        const real_t Omega = (real_t)10.0;
        const Vector3r position = Vector3r::Zero();
        const Quaternion4r orientation = Quaternion4r::Identity();
    PhysicsSystem::CubeShape shape{Vector3r(a / 2, b / 2, c / 2)};
        const Vector3r linearVelocity = Vector3r::Zero();
        // const Vector3r angularVelocity(Omega, eps, eps);
        const Vector3r angularVelocity(eps, Omega, eps);
        // const Vector3r angularVelocity(eps, eps, Omega);
    PhysicsSystem::RigidState state; state.position = position; state.orientation = orientation; state.linearVelocity = linearVelocity; state.angularVelocity = angularVelocity; PhysicsSystem::RigidProps props; props.mass = mass; sys.addRigidBody(shape, state, props);

        const auto& reg = sys.ecs();
        entt::entity rodEntity = entt::null;
        for (auto entity : reg.view<cardillo::PhysicsSystem::C_RigidBodyTag>()) {
            rodEntity = entity; break;
        }
        
        // Print Inertia and other stats from physics system:
        const real_t* inertiaDiag = reg.get<cardillo::PhysicsSystem::C_InertiaDiag>(rodEntity).I.data();
        PetscPrintf(PETSC_COMM_WORLD, "Dzhanibekov inertia diag: A=%.6f, B=%.6f, C=%.6f\n",
                    inertiaDiag[0], inertiaDiag[1], inertiaDiag[2]);
        PetscPrintf(PETSC_COMM_WORLD, " A = %.6f\n", (real_t)(1.0/12.0) * (b * b + c * c));
        PetscPrintf(PETSC_COMM_WORLD, " B = %.6f\n", (real_t)(1.0/12.0) * (a * a + c * c));
        PetscPrintf(PETSC_COMM_WORLD, " C = %.6f\n", (real_t)(1.0/12.0) * (a * a + b * b));

    }
};
