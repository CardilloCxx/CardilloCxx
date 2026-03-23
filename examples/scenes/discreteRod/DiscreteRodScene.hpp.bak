#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace cardillo;

// Interlocking half-circles scene with attached rods
class DiscreteRodScene : public SceneBase {
public:
    const char* sceneName() const override { return "discrete_rod"; }
    DiscreteRodScene() = default;
    ~DiscreteRodScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        auto& sys = engine.world();
        using namespace cardillo;
        using namespace cardillo::misc;
        // Simple ground
        engine.addStaticBody(physics::CubeShape(Vector3r(15.0, 15.0, 0.1)), physics::RigidState(Vector3r(0,0,-0.1)));

        // Common material/geometry
        const size_t seg_arc = 160;         // smooth semicircles
        const real_t radius = (real_t)0.8;  // D radius
        const real_t density = (real_t)600;
        const real_t E = (real_t)5e7;
        const real_t nu = (real_t)0.3;
        physics::BeamCrossSection section((real_t)0.01, (real_t)0.01, physics::BeamBodyType::Capsule);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);

        // Center height and layout
        const Vector3r C(0, 0, (real_t)1.5);

        // Create two rigid rods (the straight sides of the Ds), very heavy
        const Vector3r halfExtents((real_t)0.03, radius * 1.1, (real_t)0.03); // length 2R along Y
        m_leftRod = engine.addRigidBody(physics::CubeShape(halfExtents),
                                     physics::RigidState(Vector3r(C.x() - radius * 0.8, 0, C.z()),
                                       Quaternion4r(Eigen::AngleAxis<real_t>((real_t)M_PI, Vector3r::UnitY()))),
                                     physics::RigidProps((real_t)1e10));
        m_rightRod = engine.addRigidBody(physics::CubeShape(halfExtents),
                                      physics::RigidState(Vector3r(C.x() + radius * 0.8, 0, C.z()),
                                        Quaternion4r(Eigen::AngleAxis<real_t>((real_t)0.5*M_PI, Vector3r::UnitX()))),
                                      physics::RigidProps((real_t)1e10));

        // Build two half-circles in the XY-plane (normal Z), opening right and left
        const real_t thetaSpan = (real_t)M_PI;
        CircleSpline arc(Vector3r::Zero(), radius, Vector3r::UnitZ(), Vector3r::UnitX(), (real_t)(0.5*M_PI), thetaSpan);

        auto arcRightEnds = engine.createBeam(arc, section, springs, physics::RigidState(Vector3r(-halfExtents.x(), 0, 0), m_leftRod, sys.ecs()),  physics::RigidProps::withDensity(density), seg_arc);
        auto arcLeftEnds  = engine.createBeam(arc, section, springs, physics::RigidState(Vector3r(-halfExtents.x(), 0, 0), m_rightRod, sys.ecs()), physics::RigidProps::withDensity(density), seg_arc);

        // Attach semicircle endpoints to corresponding rods to form two Ds
        sys.addRigidConstraint(m_leftRod,  arcRightEnds.first);
        sys.disableCollisionBetween(m_leftRod,  arcRightEnds.first);
        sys.addRigidConstraint(m_leftRod,  arcRightEnds.second);
        sys.disableCollisionBetween(m_leftRod,  arcRightEnds.second);

        sys.addRigidConstraint(m_rightRod, arcLeftEnds.first);
        sys.disableCollisionBetween(m_rightRod, arcLeftEnds.first);
        sys.addRigidConstraint(m_rightRod, arcLeftEnds.second);
        sys.disableCollisionBetween(m_rightRod, arcLeftEnds.second);
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        auto& sys = engine.world();
        (void)t;
        // Cancel gravity and prescribe motion for the two rods
        const real_t v_sep = (real_t)0.5;    // move apart along ±X
        const real_t omega = (real_t)1.0;    // twist around local Y

        if (m_leftRod != entt::null) {
            sys.applyForce(m_leftRod, -sys.gravity() * sys.getMass(m_leftRod).diagonal(), Vector3r::Zero());
            sys.setLinearVelocity(m_leftRod, Vector3r(-(v_sep), 0, 0));
            sys.setAngularVelocity(m_leftRod, Vector3r(omega, 0, 0));
        }
        if (m_rightRod != entt::null) {
            sys.applyForce(m_rightRod, -sys.gravity() * sys.getMass(m_rightRod).diagonal(), Vector3r::Zero());
            sys.setLinearVelocity(m_rightRod, Vector3r(+(v_sep), 0, 0));
            sys.setAngularVelocity(m_rightRod, Vector3r(omega, 0, 0));
        }
    }

    private:
        entt::entity m_leftRod{entt::null};
        entt::entity m_rightRod{entt::null};
};