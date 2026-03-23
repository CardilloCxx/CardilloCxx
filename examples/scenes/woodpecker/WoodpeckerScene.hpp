#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <iostream>
#include <cmath>

// Woodpecker scene: a static vertical rod, a ring mesh, and a
// woodpecker mesh connected via a hinge at the ring's edge.
class WoodpeckerScene : public SceneBase {
public:
    const char* sceneName() const override { return "woodpecker"; }
    WoodpeckerScene() = default;
    ~WoodpeckerScene() override = default; 

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        std::cout << "Populating Woodpecker scene..." << std::endl;

        // Toggle: keep the original compliant joint constraint instead of the helical spring.
        // (Set to true to restore the previous behavior.)
        const bool keepOldConstraint = false;

        // 1) Static vertical rod as a capsule: radius 1.7mm, height 30cm (halfLength = 0.15)
        const real_t rodRadius = (real_t)0.0017;
        const real_t rodHalfLength = (real_t)0.15; // 30cm total
        physics::CapsuleShape rodShape(rodRadius, rodHalfLength);
        physics::RigidState rodState(Vector3r(0.0,0.0,-0.145));
        engine.addStaticBody(rodShape, rodState);

        // 2) Dynamic ring and woodpecker meshes (already authored at correct positions)
        //    Use modest densities to compute mass from volume automatically
        const real_t rhoRing = (real_t)400;      // kg/m^3
        const real_t rhoWoodpecker = (real_t)400; // kg/m^3
        auto ring = engine.addRigidBody(physics::MeshShape("res/meshes/woodpecker_ring.obj"),
                                     physics::RigidState(Vector3r::Zero(), Vector3r(0.0, 0.0, 0.0)),
                                     physics::RigidProps::withDensity(rhoRing));

        auto pecker = engine.addRigidBody(physics::MeshShape("res/meshes/woodpecker.obj"),
                                       physics::RigidState(Vector3r::Zero(), Vector3r(0.0, 0.0, -0.2)),
                                       physics::RigidProps::withDensity(rhoWoodpecker));

        // 3) Place hinge at ring edge (radius 5mm) in -Y direction and offset woodpecker 1cm further -Y
        const real_t ringEdgeRadius = (real_t)0.005; // 5mm
        const real_t woodpeckerOffset = (real_t)0.01; // 1cm from ring edge towards -Y

        // 4) Connection between ring and woodpecker:
        //    Either keep the old compliant joint (TranslationRotationConstraint)
        //    or replace it by a helical beam spring with fixed rest-length and fixed turns.
        const Vector3r hingeCenter = Vector3r(0, -ringEdgeRadius, 0);
        if (keepOldConstraint) {
            physics::JointFrame jf = physics::JointFrame(hingeCenter, Matrix33r::Identity(), std::nullopt);

            engine.addTranslationRotationConstraint(ring, pecker, jf,
                                               /*K_trans*/ Vector3r(3000, std::numeric_limits<real_t>::infinity(), std::numeric_limits<real_t>::infinity()),
                                               /*D_trans*/ Vector3r::Zero(),
                                               /*K_rot*/   Vector3r(0.05, 0.5, 0.05),
                                               /*D_rot*/   Vector3r::Zero());
        } else {
            using namespace cardillo;
            using namespace cardillo::misc;

            const real_t wireDiameter = (real_t)0.0003;                // 0.3 mm
            const real_t coilRadius   = (real_t)0.001;                  // 1 mm
            const size_t turns        = 8;                             // fixed turns
            const size_t segmentsPerTurn = 24;
            const size_t segments = turns * segmentsPerTurn;
            const Vector3r pA = hingeCenter + Vector3r(0, 0, coilRadius * 0.5);
            const Vector3r pB = hingeCenter + Vector3r(0, -woodpeckerOffset, coilRadius * 0.5);
            const real_t density = (real_t)7850.0;
            const real_t E = (real_t)206e9;
            const real_t nu = (real_t)0.2638;

            const Vector3r axis = (pB - pA);
            const real_t axisLen = axis.norm();
            Vector3r axisDir = Vector3r::UnitX();
            if (axisLen > (real_t)1e-12) axisDir = axis / axisLen;

            const real_t pitch = axisLen / static_cast<real_t>(turns);

            Vector3r dir0 = Vector3r::UnitZ();
            Vector3r u = dir0 - dir0.dot(axisDir) * axisDir;
            if (u.squaredNorm() < (real_t)1e-12) {
                dir0 = Vector3r::UnitY();
                u = dir0 - dir0.dot(axisDir) * axisDir;
            }
            u.normalize();

            const Vector3r helixCenter = pA - coilRadius * u;
            auto helix = HelixSpline(helixCenter, axisDir, coilRadius, pitch, static_cast<real_t>(turns), dir0);
            std::vector<const misc::SplinePattern*> parts{&helix};

        
            physics::BeamCrossSection sec(wireDiameter, wireDiameter, physics::BeamBodyType::Capsule);
            auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
            springs.setDampingFromFactor((real_t)0.001);
            auto props = physics::RigidProps::withDensity(density);
            props.collidable = false; 

            auto endpoints = engine.createBeams(parts, sec, springs,
                                             physics::RigidState{},
                                             props,
                                             segments);

            const Vector3r inf3 = Vector3r::Constant(std::numeric_limits<real_t>::infinity());
            engine.addTranslationRotationConstraint(endpoints.first, ring,
                physics::JointFrame(pA, Matrix33r::Identity(), std::nullopt),
                /*K_trans*/ inf3, /*D_trans*/ Vector3r::Zero(),
                /*K_rot*/   inf3, /*D_rot*/   Vector3r::Zero());
            engine.addTranslationRotationConstraint(endpoints.second, pecker,
                physics::JointFrame(pB, Matrix33r::Identity(), std::nullopt),
                /*K_trans*/ inf3, /*D_trans*/ Vector3r::Zero(),
                /*K_rot*/   inf3, /*D_rot*/   Vector3r::Zero());
        }

        // 5) Track entities for output
        engine.track(pecker, "woodpecker_body");
        engine.track(ring, "woodpecker_ring");
    }
};
