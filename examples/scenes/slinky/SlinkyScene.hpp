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

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;
        // Build staircase as static cubes
        const int steps = 12;
        const real_t rise = (real_t)0.0025 * 50.0;
        const real_t run  = (real_t)0.0025 * 50.0;
        const real_t width = (real_t) 1.0;           // wider stairs
        const real_t thickness = (real_t)0.04; 
        const Vector3r stairOrigin(0, 0, 0);

        entt::entity top_step = entt::null;
        for (int i = 0; i < steps; ++i) {
            real_t z = stairOrigin.z() + (real_t)i * rise * 1.2;
            real_t x = stairOrigin.x() + (real_t)i * run;
            Vector3r halfExtents(run * (real_t)0.5, width * (real_t)0.5, thickness * (real_t)0.5);
            Vector3r pos(x + halfExtents.x(), 0, z - thickness * (real_t)0.5);
            top_step = engine.addStaticBody(physics::CubeShape(halfExtents), physics::RigidState(pos));
            // engine.addStaticBody(physics::SphereShape(thickness * (real_t)0.9), physics::RigidState(pos + Vector3r(0, 0, halfExtents.z())));
            engine.addStaticBody(physics::ConeShape(thickness * 0.9, thickness * 1.5), physics::RigidState(pos + Vector3r(0, 0, thickness * 0.75 + halfExtents.z())));
        }

        // Slinky parameters (plastic-like)
        const int turns = 50;                 // number of windings
        const real_t diameter = (real_t)0.065; // ~6.5 cm OD
        const real_t pitch = (real_t)rise / turns;   // axial rise per turn ~2.5 mm
        const real_t radius = diameter * (real_t)0.7;
        const size_t segments = turns * 16.5;   // uneven per turn to avoid interlocking
        const real_t density = (real_t)700;   // plastic density kg/m^3
        const real_t E = (real_t)1e9;         // lower Young's modulus for plastic 1e9, 5e9 for metal
        const real_t nu = (real_t)0.35;
        physics::BeamCrossSection section(pitch*2, pitch*0.99, physics::BeamBodyType::Cube); // (pitch*0.99, pitch*0.99, physics::BeamBodyType::Capsule); for metal pitch*2, pitch*0.99, physics::BeamBodyType::Cube for plastic-like
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        springs.setDampingFromFactor(0.00);

        // Place the spring centered over the top step tread and build a helix
        const real_t topZ = stairOrigin.z() + (real_t)(steps - 1) * rise;
        const real_t topX = stairOrigin.x() + (real_t)(steps - 1) * run + run * (real_t)0.5; // tread center
        Vector3r topStepCenter(topX, (real_t)0.0, topZ);
        const real_t lift = (real_t)0.01; // small clearance above tread
        Vector3r center0 = topStepCenter + Vector3r((real_t)0.0, (real_t)0.0, lift);

        misc::HelixSpline helix(Vector3r::Zero(), Vector3r::UnitZ(), radius, pitch, (real_t)turns, Vector3r::UnitX());

        std::cout << "Creating slinky with " << segments << " segments, " << turns << " turns, pitch " << pitch << ", radius " << radius << "\n";
        physics::RigidState stateDefaults(Vector3r(0.0,0.0,thickness * 0.5 + pitch * 0.5), Vector3r::Zero(),  Quaternion4r::Identity(), Vector3r(0.0,0.0,0.0), top_step, engine.ecs());
        physics::RigidProps props = physics::RigidProps::withDensity(density);
        auto ends = engine.createBeam(helix, section, springs, stateDefaults, props, segments);
        m_slinky_end_entity = ends.second;
    
        // m_guide_entity = engine.addRigidBody(physics::SphereShape(radius * (real_t)0.001), 
        //                                   physics::RigidState(ends.second), 
        //                                   physics::RigidProps(1e10));

        // engine.setLinearVelocity(ends.second, Vector3r((real_t)-5.0, (real_t)0.0, (real_t)5.0));
        // engine.addRigidConstraint(m_guide_entity,  ends.second);
        // engine.disableCollisionBetween(m_guide_entity, ends.first);
    
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t dt) override {
        (void)dt;
        auto A_IK = engine.ecs().get<C_Orientation>(m_slinky_end_entity).value.toRotationMatrix();
        
        if (t < 0.28) engine.applyForce(m_slinky_end_entity, Vector3r(0.00, 0.0, 0.0), A_IK.transpose()* Vector3r(0.0, -0.015, 0.0));
//         engine.setAngularVelocity(m_guide_entity, Vector3r(0.0, -0.1, 0.0));
//         engine.applyForce(m_guide_entity, -engine.gravity() * engine.getMass(m_guide_entity).diagonal(), Vector3r::Zero());
    }

private: 
    entt::entity m_slinky_end_entity;
    entt::entity m_guide_entity;
};
