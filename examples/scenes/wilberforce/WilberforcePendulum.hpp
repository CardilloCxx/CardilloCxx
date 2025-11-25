#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <memory>
#include <cmath>

using namespace cardillo;

// Wilberforce Pendulum:
// A mass-spring system where torsional and vertical modes are tuned to exchange energy.
// We model the spring as a helical beam (capsule cross-section) anchored at the top,
// with a rigid mass at the bottom. Initial conditions excite vertical mode and small twist.
class WilberforcePendulumScene : public SceneBase {
public:
    const char* sceneName() const override { return "wilberforce"; }
    WilberforcePendulumScene() = default;
    ~WilberforcePendulumScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Spring parameters
        const real_t wireRadius = 0.001;              // wire thickness    1  mm
        const real_t meanRadius = 0.015;              // spring mean radius 15 mm
        const size_t turns = 50;                      // helical turns count
        const size_t segments = 1000;                 // discretization along the helix
        const real_t freeLength = wireRadius * turns * 3; 

        // Material properties
        const real_t density = 7800.0;               // kg/m^3
        const real_t E = 2.1e11;                     // Young's modulus
        const real_t nu = 0.29;                      // Poisson ratio
        const real_t G = E / (2.0 * (1.0 + nu));     // Shear modulus

        // Beam cross-section (capsule) used by createBeam
        PhysicsSystem::BeamCrossSection sec(wireRadius*2, wireRadius*2, PhysicsSystem::BeamBodyType::Capsule);
        auto springs = PhysicsSystem::BeamSpringParams::fromMaterial(E, nu);
        const real_t pitch = freeLength / static_cast<real_t>(turns);
        // Build spring via CompoundSpline (single helix segment now, extensible for future pieces)
        auto helix = misc::HelixSpline(Vector3r::Zero(), -Vector3r::UnitZ(), meanRadius, pitch, static_cast<real_t>(turns), Vector3r::UnitX());
        auto linear = misc::LinearSpline(Vector3r(meanRadius, 0, -freeLength), Vector3r(0, 0, -freeLength));
        auto linear2 = misc::LinearSpline( Vector3r(0, 0, -freeLength), Vector3r(0, 0, -1.2 * freeLength));

        // Build sequence of splines; create beams per spline and connect with rigid constraints.
        std::vector<const misc::SplinePattern*> parts{&helix, &linear, &linear2};
        auto endpoints = sys.createBeams(parts, sec, springs, PhysicsSystem::RigidState{}, PhysicsSystem::RigidProps::withDensity(density), segments);
        m_top = endpoints.first; sys.makeStatic(m_top);
        m_bottom = endpoints.second;

        // we want K / m = lambda / Iz
        const real_t d = wireRadius * 2.0;  // wire diameter
        const real_t D = meanRadius * 2.0;  // mean coil diameter
        const real_t n = static_cast<real_t>(turns); // number of active coils
        const real_t K =      (G * std::pow(d, 4)) / (8.0 * std::pow(D, 3) * n); // helix axial stiffness
        const real_t lambda = (G * std::pow(d, 4) )/ (32.0 * n * D);        // helix torsional stiffness
        const real_t tunedMass = 0.15;
        const real_t Iz = lambda / (K / tunedMass); 
        const real_t tunedSize = std::sqrt(6 * Iz / tunedMass); // Iz = 1 / 12 m (a^2 + a^2) = 1/6 m a^2 => a = sqrt(6 Iz / m)

        std::cout << "Tuning info: K = " << K << " N/m, lambda = " << lambda << " Nm/rad, m = " << tunedMass
                  << " kg, Iz = " << Iz << " kg m^2, cube size = " << tunedSize << " m" << "Expected frequencies: f_z = " << std::sqrt(K / tunedMass) / (2 * M_PI) << " Hz"
                  << " f_theta = " << std::sqrt(lambda / Iz) / (2 * M_PI) << " Hz, Spring mass: " << sys.getMass(m_bottom).col(0).row(0) * segments << " kg" << std::endl;

        m_bob = sys.addRigidBody(PhysicsSystem::CubeShape(Vector3r(tunedSize,tunedSize,tunedSize)),
        PhysicsSystem::RigidState(Vector3r(0,0,-tunedSize) + sys.getPosition(m_bottom).head<3>()),
        PhysicsSystem::RigidProps(tunedMass));

        // Pin bottom endpoint to bob using a rigid constraint
        sys.addConstraint<physics::RigidConstraint>(sys.ecs(), m_bottom, m_bob, Vector3r::Zero(), Vector3r(0,0,tunedSize));
        
        const real_t vz0 = -1.0; // m/s downward
        const real_t wz0 = 0.0;  // rad/s small spin around Z to couple torsion & vertical modes

        sys.setLinearVelocity(m_bob, Vector3r(0,0,vz0));
        sys.setAngularVelocity(m_bob, Vector3r(0,0,wz0));
    }

    void updateScene(cardillo::PhysicsSystem& /*sys*/, real_t /*t*/, real_t /*dt*/) override {

    }

private:
    entt::entity m_top{entt::null};
    entt::entity m_bottom{entt::null};
    entt::entity m_bob{entt::null};
};
