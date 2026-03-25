#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Dzhanibekov rod scene (moved into its own folder)
class DzhanibekovScene : public SceneBase {
public:
    const char* sceneName() const override { return "dzhanibekov"; }
    DzhanibekovScene() = default;
    ~DzhanibekovScene() override = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        const real_t mass = (real_t)1.0;
        const real_t a = (real_t)3.0;
        const real_t b = (real_t)2.0;
        const real_t c = (real_t)1.0;
        const real_t eps = (real_t)1e-3;
        const real_t Omega = (real_t)10.0;
        const Vector3r position = Vector3r::Zero();
        const Quaternion4r orientation = Quaternion4r::Identity();
    physics::CubeShape shape{Vector3r(a / 2, b / 2, c / 2)};
        const Vector3r linearVelocity = Vector3r::Zero();
        // const Vector3r angularVelocity(Omega, eps, eps);
        const Vector3r angularVelocity(eps, Omega, eps);
        // const Vector3r angularVelocity(eps, eps, Omega);
    physics::RigidState state; state.position = position; state.orientation = orientation; state.linearVelocity = linearVelocity; state.angularVelocity = angularVelocity; physics::RigidProps props; props.mass = mass; engine.addRigidBody(shape, state, props);

        const auto& reg = engine.ecs();
        entt::entity rodEntity = entt::null;
        for (auto entity : reg.view<cardillo::C_RigidBodyTag>()) {
            rodEntity = entity; break;
        }
        if (rodEntity != entt::null) {
            engine.track(rodEntity, "dzhanibekov_rod");
        }

    }
};
