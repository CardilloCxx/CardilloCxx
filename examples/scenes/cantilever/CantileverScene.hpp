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

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        engine.setGravity(Vector3r(0, 0, -9.81)); // no gravity for this scene

        // Material and geometry
        const real_t L  = (real_t)0.31415; // 31 cm
        const real_t r  = (real_t)0.01;    // 10 mm radius
        const real_t d  = (real_t)2 * r;   // 20 mm diameter
        const real_t E  = (real_t)7e6;     // 7 GPa
        const real_t nu = (real_t)0.75;    // Poisson ratio
        // const real_t rho = (real_t)0.4 / (M_PI * r * r * L);  // kg/m^3 (mass 0.4 kg)
        const real_t rho = (real_t)1350.0;  // kg/m^3 (typical plastic)

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

        const auto section = physics::BeamCrossSection::triangle(d,d);
        auto beam_ends = engine.createBeam(spline, section, springs, stateDefaults, props, segments, physics::BeamColliderMode::InterSegmentHull);
        m_beamRightEnd1 = beam_ends.second;
        engine.makeStatic(beam_ends.first);
        engine.addTranslationalConstraint(beam_ends.second, entt::null, physics::JointFrame(beam_ends.second));


        // Create a second beam
        const auto section2 = physics::BeamCrossSection::square(0.75* d, 0.1 * d);
        stateDefaults.setOrientation(Quaternion4r(Eigen::AngleAxis<real_t>(M_PI_2 / 2.0, Vector3r::UnitZ())));

        stateDefaults.angularVelocity = Vector3r(0, 0, 10.0);
        misc::LinearSpline spline2(Vector3r(L * 0.5, -L * 0.5, 0.1), Vector3r(L * 0.5, L * 0.5, 0.1));
        auto springs2 = physics::BeamSpringParams::fromMaterial(E * 0.4, nu);
        beam_ends = engine.createBeam(spline2, section2, springs2, stateDefaults, props, segments * 2, physics::BeamColliderMode::InterSegmentHull);

        // spheres falling on the beam
        // engine.addRigidBody(physics::SphereShape(0.01), physics::RigidState(Vector3r(L * 0.5, 0.1, 0.06)), physics::RigidProps::withDensity(20000.0));

        m_Kf = springs.Kf(L / segments, section);
        m_L = L;
        m_segments = segments;

    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        real_t t1 = 3.0;
        if (t < t1 * 10.0) {
            // engine.applyForce(m_beamRightEnd, Vector3r(-0.5 * std::max(t1, t), -0.5 * std::max(t1, t), 1.5 * std::max(t1, t)), Vector3r::Zero());
            
            real_t time = std::min(t, t1);
            Vector3r force(0.0, -1 * time / t1, 0.0);
            // Vector3r moment(0.0, -2 * M_PI * m_Kf(1) * m_L / m_segments / m_L * t / t1 / 2, 0.0);
            Vector3r torsion(5.0 * m_Kf(2) * m_L / m_segments / m_L * time / t1, 0.0, 0.0);
            engine.applyForce(m_beamRightEnd1, force * 0.0, torsion);

            // engine.applyForce(m_beamRightEnd3, force, moment);
            // engine.applyForce(m_beamRightEnd, Vector3r(0, -0.5 * std::max(t1, t), 0), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(-1.0 * std::max(t1, t), 0, 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, -1.0 * std::max(t1, t), 0));
            // engine.applyForce(m_beamRightEnd, Vector3r::Zero(), Vector3r(0, 0, -1.0 * std::max(t1, t)));
        }

        // if (t > 0.1) 
        // {
        //     // Make every entity static except the sphere:
        //     auto rb_view = engine.ecs().view<C_PhysicsObject>();
        //     for (auto [entity] : rb_view.each()) {
        //         if (entity == m_beamRightEnd1 || entity == m_beamRightEnd2 || entity == m_beamRightEnd3) continue;
        //         engine.makeStatic(entity);
        //     }
        // }
    }

    private:
     entt::entity m_beamRightEnd1{entt::null};
     Vector3r m_Kf{Vector3r::Zero()};
     real_t m_L{0.0};
     size_t m_segments{0};
};