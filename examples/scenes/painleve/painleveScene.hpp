#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Painlevé rod scene (moved into its own folder)
class PainleveScene : public SceneBase {
public:
    const char* sceneName() const override { return "painleve"; }
    PainleveScene() = default;
    ~PainleveScene() override = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        auto& sys = engine.world();
        using namespace cardillo;

        const real_t groundHalfThickness = (real_t)0.01;
        const real_t groundHalfSize = (real_t)6000.0;
        physics::CubeShape groundShape{Vector3r(groundHalfSize, groundHalfSize, groundHalfThickness)};
        physics::RigidState groundState; groundState.position = Vector3r(0.0,0.0,-groundHalfThickness); groundState.orientation = Quaternion4r::Identity();
        engine.addStaticBody(groundShape, groundState);

        const real_t phi = (real_t)(-31.0 * M_PI / 180.0);
        const real_t friction = (real_t)5 / 2;
        const real_t mass = (real_t)1.0;
        const real_t J = (real_t)(mass * (1.0 / 3.0));
        const real_t l = (real_t)1.0;
        const real_t thickness = (real_t)0.02;

        // const real_t R = thickness * (real_t)0.5;
        // const real_t z_com = l * std::sin(-phi) + R; 
        // const Vector3r position(0.0, 0.0, z_com);
        // physics::CapsuleShape rodShape{thickness * (real_t)0.5, l};
        // const Quaternion4r capsuleOrientation(Eigen::AngleAxis<real_t>(M_PI_2 -phi, Vector3r::UnitY()));
        // const real_t v0 = (real_t)30.0;
        // const Vector3r linearVelocity(v0, 0.0, 0.0);
        // physics::RigidState state; state.position = position; state.orientation = capsuleOrientation; state.linearVelocity = linearVelocity; state.angularVelocity = Vector3r::Zero();
        // physics::RigidProps props; props.mass = mass; props.friction = friction;
        // auto rod = engine.addRigidBody(rodShape, state, props);
        // sys.ecs().get<World::C_InertiaDiag>(rod).I = Vector3r(J, J, J);
        // sys.track(rod, "painleve_rod");

        // Parameter sweep over friction and initial angle phi
        const size_t num_friction_steps = 100;
        const real_t friction_min = 0.5;                    const real_t friction_max = 5.5;     
         
        const size_t num_phi_steps = 100;
        const real_t phi_min =  (real_t)(10 * M_PI / 180.0); const real_t phi_max = (real_t)(90 * M_PI / 180.0);

        const real_t spacing = 5.0;

        for (size_t i = 0; i < num_friction_steps; ++i) {
            real_t curr_friction = friction_min + (friction_max - friction_min) * ((real_t)i / (real_t)(num_friction_steps - 1));
            for (size_t j = 0; j < num_phi_steps; ++j) {
                real_t curr_phi = phi_min + (phi_max - phi_min) * ((real_t)j / (real_t)(num_phi_steps - 1));

                const real_t x = (real_t)spacing * (i );
                const real_t y = (real_t)spacing * (j );
                const real_t R = thickness * (real_t)0.5;
                const real_t z_com = l * std::sin(curr_phi) + R; 
                const Vector3r position(x, y, z_com);
                physics::CapsuleShape rodShape{thickness * (real_t)0.5, l};
                const Quaternion4r capsuleOrientation(Eigen::AngleAxis<real_t>(M_PI_2 -curr_phi, Vector3r::UnitY()));
                const real_t v0 = (real_t)30.0;
                const Vector3r linearVelocity(-v0, 0.0, 0.0);
                physics::RigidState state; state.position = position; state.orientation = capsuleOrientation; state.linearVelocity = linearVelocity; state.angularVelocity = Vector3r::Zero();
                physics::RigidProps props; props.mass = mass; props.friction = curr_friction;
                auto rod = engine.addRigidBody(rodShape, state, props);
                sys.ecs().get<World::C_InertiaDiag>(rod).I = Vector3r(J, J, J);
                sys.track(rod, "painleve_rod_fric_" + std::to_string(curr_friction) + "_phi_" + std::to_string(curr_phi));
            }
        }
    }
};
