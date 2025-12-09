#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

using namespace cardillo;

// Single spaghetti modeled as a straight beam made of capsule segments
class SpaghettiScene : public SceneBase {
public:
    const char* sceneName() const override { return "spaghetti"; }
    SpaghettiScene() = default;
    ~SpaghettiScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;
    
        sys.setGravity(Vector3r(0,0, 9.81)); // no gravity for this scene

        // Material and geometry for spaghetti
        const real_t E  = (real_t)3.8e9;       // 3.8 GPa
        const real_t nu = (real_t)0.3;         // Poisson ratio
        const real_t rho = (real_t)1500.0;     // 1.5 g/cm^3 = 1500 kg/m^3
        const real_t L  = (real_t)0.240;       // 240 mm
        const real_t d  = (real_t)0.00075;     // 0.75 mm

        // Cross-section: circular (capsule body type uses circle properties via min(width,height))
        PhysicsSystem::BeamCrossSection section(d, d, PhysicsSystem::BeamBodyType::Capsule);

        // Beam spring params from material (no extra damping by default).
        auto springs = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);
        springs.tensileStrength = 20e6; // 50 MPa tensile strength
        springs.crackStrainMax = 0.001; // shear retention threshold

        // Build a straight beam along +X of length L, elevated above ground
        const Vector3r p0(0, 0, (real_t)0.05 + d);
        const Vector3r p1 = p0 + Vector3r(L, 0, 0);
        misc::LinearSpline spline(p0, p1);

        // Choose a reasonable segment count so segment length ~ few diameters for stability
        const real_t segLenTarget = d * (real_t)1.0; // ~6 diameters
        const size_t segments = (size_t)std::max<real_t>(1000, std::ceil(L / segLenTarget));

        // Default state: identity orientation; default density props
        PhysicsSystem::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());
        PhysicsSystem::RigidProps props = PhysicsSystem::RigidProps::withDensity(rho);

        auto beam_ends = sys.createBeam(spline, section, springs, stateDefaults, props, segments);
  

        m_beamLeftEnd = sys.addRigidBody(PhysicsSystem::SphereShape(d * (real_t)0.0001),
                          PhysicsSystem::RigidState(p0 - Vector3r(d,0,0)),
                          PhysicsSystem::RigidProps(1e10));
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_beamLeftEnd, beam_ends.first, physics::JointFrame(m_beamLeftEnd));

        m_beamRightEnd = sys.addRigidBody(PhysicsSystem::SphereShape(d * (real_t)0.0001),
                           PhysicsSystem::RigidState(p1 + Vector3r(d,0,0)),
                           PhysicsSystem::RigidProps(1e10));
        sys.addConstraint<physics::TranslationalConstraint>(sys.ecs(), m_beamRightEnd, beam_ends.second, physics::JointFrame(m_beamRightEnd));

    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t /*t*/, real_t /*dt*/) override {
        
        const real_t speed = 0.02; // m/s
        sys.setLinearVelocity(m_beamLeftEnd, Vector3r(speed,0,0));
        sys.setLinearVelocity(m_beamRightEnd, Vector3r(-speed,0,0));

        sys.applyForce(m_beamLeftEnd, 1e10 * 9.81 * Vector3r(0,0,-1), Vector3r::Zero());
        sys.applyForce(m_beamRightEnd, 1e10 * 9.81 * Vector3r(0,0,-1), Vector3r::Zero());
    }

private:
    entt::entity m_beamLeftEnd{entt::null};
    entt::entity m_beamRightEnd{entt::null};
};
