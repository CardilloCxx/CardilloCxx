#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <vector>
#include <cmath>

using namespace cardillo;

class SlinkyScene : public SceneBase {
public:
    const char* sceneName() const override { return "slinky"; }
    SlinkyScene() = default;
    ~SlinkyScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;
        // Build staircase as static cubes
        const int steps = 12;
        const real_t rise = (real_t)0.0025 * 50.0;
        const real_t run  = (real_t)0.0025 * 50.0;
        const real_t width = (real_t) 1.0;           // wider stairs
        const real_t thickness = (real_t)0.04; 
        const Vector3r stairOrigin(0, 0, 0);

        entt::entity top_step = entt::null;
        for (int i = 0; i < steps; ++i) {
            real_t z = stairOrigin.z() + (real_t)i * rise;
            real_t x = stairOrigin.x() + (real_t)i * run;
            Vector3r halfExtents(run * (real_t)0.5, width * (real_t)0.5, thickness * (real_t)0.5);
            Vector3r pos(x + halfExtents.x(), 0, z - thickness * (real_t)0.5);
            top_step = sys.addStaticBody(PhysicsSystem::CubeShape(halfExtents), PhysicsSystem::RigidState(pos));
            sys.addStaticBody(PhysicsSystem::SphereShape(thickness * (real_t)0.9), PhysicsSystem::RigidState(pos + Vector3r(0, 0, halfExtents.z())));
        }

        // Slinky parameters (plastic-like)
        const int turns = 50;                 // number of windings
        const real_t diameter = (real_t)0.065; // ~6.5 cm OD
        const real_t pitch = (real_t)rise / turns;   // axial rise per turn ~2.5 mm
        const real_t radius = diameter * (real_t)0.7;
        const size_t segments = turns * 16.5;   // fewer segments per turn for performance
        const real_t density = (real_t)700;   // plastic density kg/m^3
        const real_t E = (real_t)5e9;         // lower Young's modulus for plastic 1e9 
        const real_t nu = (real_t)0.35;
        PhysicsSystem::BeamCrossSection section(pitch*0.99, pitch*0.99, PhysicsSystem::BeamBodyType::Capsule); // pitch*2, pitch*0.99, PhysicsSystem::BeamBodyType::Cube for plastic-like
        auto springs = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);

        // Place the spring centered over the top step tread and build a helix
        const real_t topZ = stairOrigin.z() + (real_t)(steps - 1) * rise;
        const real_t topX = stairOrigin.x() + (real_t)(steps - 1) * run + run * (real_t)0.5; // tread center
        Vector3r topStepCenter(topX, (real_t)0.0, topZ);
        const real_t lift = (real_t)0.01; // small clearance above tread
        Vector3r center0 = topStepCenter + Vector3r((real_t)0.0, (real_t)0.0, lift);

        misc::HelixSpline helix(Vector3r::Zero(), Vector3r::UnitZ(), radius, pitch, (real_t)turns, Vector3r::UnitX());

        std::cout << "Creating slinky with " << segments << " segments, " << turns << " turns, pitch " << pitch << ", radius " << radius << "\n";
        PhysicsSystem::RigidState stateDefaults(Vector3r(0.0,0.0,thickness * 0.5 + pitch * 0.5), Vector3r::Zero(),  Quaternion4r::Identity(), Vector3r(0.0,0.0,0.0), top_step, sys.ecs());
        PhysicsSystem::RigidProps props = PhysicsSystem::RigidProps::withDensity(density);
        auto ends = sys.createBeam(helix, section, springs, stateDefaults, props, segments);
        m_slinky_end_entity = ends.second;
    
        // m_guide_entity = sys.addRigidBody(PhysicsSystem::SphereShape(radius * (real_t)0.001), 
        //                                   PhysicsSystem::RigidState(ends.second), 
        //                                   PhysicsSystem::RigidProps(1e10));

        // sys.setLinearVelocity(ends.second, Vector3r((real_t)-5.0, (real_t)0.0, (real_t)5.0));
        // sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_guide_entity,  ends.second);
        // sys.disableCollisionBetween(m_guide_entity, ends.first);
    
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t dt) override {
        auto A_IK = sys.ecs().get<cardillo::PhysicsSystem::C_Orientation>(m_slinky_end_entity).value.toRotationMatrix();
        
        if (t < 0.15) sys.applyForce(m_slinky_end_entity, Vector3r(-0.1, 0.0, 0.0), A_IK.transpose()* Vector3r(0.0, -0.02, 0.0));
//         sys.setAngularVelocity(m_guide_entity, Vector3r(0.0, -0.1, 0.0));
//         sys.applyForce(m_guide_entity, -sys.gravity() * sys.getMass(m_guide_entity).diagonal(), Vector3r::Zero());
    }

private: 
    entt::entity m_slinky_end_entity;
    entt::entity m_guide_entity;
};
