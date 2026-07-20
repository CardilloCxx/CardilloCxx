#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

class NetScene : public SceneBase {
public:
    const char* sceneName() const override { return "net"; }
    NetScene() = default;
    ~NetScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        const real_t radius = 0.5;
        const real_t length = 2.0 * M_PI * radius;
        const real_t thickness = 0.03;

        // Steel properties
        const real_t density = 7850; // kg/m^3
        const real_t E = 2.1e11;
        const real_t nu = 0.3;

        misc::CircleSpline ring(Vector3r::Zero(), radius, Vector3r::UnitX(), Vector3r::UnitZ());
        physics::BeamCrossSection sec_ring(thickness, thickness, physics::BeamBodyType::Capsule); 
        auto springs_ring = physics::BeamSpringParams::fromMaterial(E, nu);

        const int n = 12;
        const int m = 18;
        const size_t segments = 32;
        const real_t eps = 1e-3;

        const real_t c = (radius * 2 - thickness * 1.5);
        const real_t a = (radius + thickness * 0.5 + eps);
        const real_t b = std::sqrt(c*c - a *a);
        const real_t pitch = (real_t) M_PI / 12.0;
        Vector3r netCenter = Vector3r(0.0, 0.0, 10.0);
        Vector3r offset = netCenter + Vector3r(0.0, -a * n, - b * m * 0.5);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                Vector3r center = offset + Vector3r(0, a * 2 * i + ((j % 2 == 0) ? 0 : a), b * j);
                const real_t angle = (j % 2 == 0) ? -pitch : pitch;
                Quaternion4r rotation = Quaternion4r(Eigen::AngleAxis<real_t>(angle, Vector3r::UnitY()));
                physics::RigidProps props;

                if ((i > 0 || j % 2 == 1) && (i < n-1 || j % 2 == 0) && (j > 0 && j < m-1)) props = physics::RigidProps::withDensity(density);
                else props = physics::RigidProps(0); // static boundary rings

                engine.createBeam(ring, sec_ring, springs_ring, physics::RigidState{center, rotation}, props, segments);
            }
        }

        // Boulder
        auto boulder = engine.addRigidBody(physics::MeshShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/rock.obj", Vector3r(0.75,0.75,0.75)),
                         physics::RigidState(Vector3r(-3,0,0) + netCenter, Vector3r(30.0, 0.0, 0.0), Vector3r(10,20,50)),
                         physics::RigidProps::withDensity(2500));

        std::cout << "Boulder mass: " << engine.getMass(boulder).col(0).row(0) << " kg" << std::endl;
        std::cout << "Boulder KE: " << engine.getKineticEnergy(boulder) << " J" << std::endl;
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        (void)engine;
        (void)t;
//         // Apply a twisting moment at the rod end
//         real_t torque_magnitude = 0.05;
//         engine.applyForce(m_rodEnd, Vector3r::Zero(), Vector3r(0, -torque_magnitude, 0));
//         // engine.applyForce(m_ropeEnd, Vector3r::Zero(), Vector3r(0.1, 0, 0));
// 
//         if (t < 0.05) {
//             engine.setLinearVelocity(m_endCube, Vector3r(0, -1, 0));
//             engine.setAngularVelocity(m_endCube, Vector3r::Zero());
//         }else {
//             engine.setLinearVelocity(m_endCube, Vector3r::Zero());
//             engine.setAngularVelocity(m_endCube, Vector3r(0, 2, 0));
//         }
//         engine.applyForce(m_endCube, -engine.gravity() * engine.getMass(m_endCube).diagonal(), Vector3r::Zero());
    }

    private:
        entt::entity m_rodEnd{entt::null};
        entt::entity m_ropeEnd{entt::null};
        entt::entity m_endCube{entt::null};
};
