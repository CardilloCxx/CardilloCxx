#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace cardillo;

// vibrations of a cantilever beam
class CantileverScene : public SceneBase {
public:
    const char* sceneName() const override { return "cantilever"; }
    CantileverScene() = default;
    ~CantileverScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        engine.setGravity(Vector3r(0, 0, 0)); // no gravity for this scene

        // Material and geometry for spaghetti
        const real_t L  = (real_t)0.31415; // 31 cm
        const real_t r  = (real_t)0.01;    // 10 mm radius
        const real_t d  = (real_t)2 * r;   // 20 mm diameter
        const real_t E  = (real_t)7e5;     // 7 GPa
        const real_t nu = (real_t)0.75;    // Poisson ratio
        // const real_t rho = (real_t)0.4 / (M_PI * r * r * L);  // kg/m^3 (mass 0.4 kg)
        const real_t rho = (real_t)1350.0;  // kg/m^3 (typical plastic)

        // Cross-section: circular (capsule body type uses circle properties via min(width,height))
        // physics::BeamCrossSection section(d, d, physics::BeamBodyType::Capsule);
        physics::BeamCrossSection section(d, d, physics::BeamBodyType::Cube);

        // Beam spring params from material (no extra damping by default).
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);

        // Build a straight beam along +X of length L, elevated above ground
        const Vector3r p0(0, 0, 0);
        const Vector3r p1 = p0 + Vector3r(L, 0, 0);
        misc::LinearSpline spline(p0, p1);

        // Number of beam segments
        const size_t segments = 20;

        // Default state: identity orientation; default density props
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());
        physics::RigidProps props = physics::RigidProps::withDensity(rho);

        auto beam_ends = engine.createBeam(spline, section, springs, stateDefaults, props, segments);
        m_beamRightEnd1 = beam_ends.second;
        engine.makeStatic(beam_ends.first);

        stateDefaults.position += 0.1 * Vector3r::UnitY();  // small offset to avoid initial self-collision
        beam_ends = engine.createBeam(spline, section, springs, stateDefaults, props, segments * 2);
        m_beamRightEnd2 = beam_ends.second;
        engine.makeStatic(beam_ends.first);

        stateDefaults.position += 0.1 * Vector3r::UnitY();  // small offset to avoid initial self-collision
        beam_ends = engine.createBeam(spline, section, springs, stateDefaults, props, segments * 4);
        m_beamRightEnd3 = beam_ends.second;
        engine.makeStatic(beam_ends.first);

        m_Kf = springs.Kf(L / segments, section);
        m_L = L;
        m_segments = segments;

        // // engine.disableCollisionBetween(m_cube, cube2);
        // auto view = engine.ecs().view<World::C_Collidable>();
        // for (auto e : view.each()) {
        //     engine.ecs().remove<World::C_PhysicsObject>(e);
        // }
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        real_t t1 = 5.0;
        if (t < t1) {
            // engine.applyForce(m_beamRightEnd, Vector3r(-0.5 * std::max(t1, t), -0.5 * std::max(t1, t), 1.5 * std::max(t1, t)), Vector3r::Zero());
            Vector3r force(0.0, -0.1 * t / t1, 0.0);
            Vector3r moment(0.0, -2 * M_PI * m_Kf(1) * m_L / m_segments / m_L * t / t1 / 2, 0.0);
            engine.applyForce(m_beamRightEnd1, force, moment);
            engine.applyForce(m_beamRightEnd2, force, moment);
            engine.applyForce(m_beamRightEnd3, force, moment);
            // engine.applyForce(m_beamRightEnd, Vector3r(0, -0.5 * std::max(t1, t), 0), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(-1.0 * std::max(t1, t), 0, 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, 0, -1.0 * std::max(t1, t)));
        }
    }

    private:
     entt::entity m_beamRightEnd1{entt::null};
     entt::entity m_beamRightEnd2{entt::null};
     entt::entity m_beamRightEnd3{entt::null};
     Vector3r m_Kf{Vector3r::Zero()};
     real_t m_L{0.0};
     size_t m_segments{0};
};