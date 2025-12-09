#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <iostream>

// Woodpecker scene: a static vertical rod, a ring mesh, and a
// woodpecker mesh connected via a hinge at the ring's edge.
class WoodpeckerScene : public SceneBase {
public:
    const char* sceneName() const override { return "woodpecker"; }
    WoodpeckerScene() = default;
    ~WoodpeckerScene() override = default; 

    void populate(cardillo::PhysicsSystem& sys) override {
        std::cout << "Populating Woodpecker scene..." << std::endl;

        // 1) Static vertical rod as a capsule: radius 1.7mm, height 30cm (halfLength = 0.15)
        const real_t rodRadius = (real_t)0.0017;
        const real_t rodHalfLength = (real_t)0.15; // 30cm total
        PhysicsSystem::CapsuleShape rodShape(rodRadius, rodHalfLength);
        PhysicsSystem::RigidState rodState(Vector3r(0.0,0.0,-0.145));
        sys.addStaticBody(rodShape, rodState);

        // 2) Dynamic ring and woodpecker meshes (already authored at correct positions)
        //    Use modest densities to compute mass from volume automatically
        const real_t rhoRing = (real_t)1500;      // kg/m^3
        const real_t rhoWoodpecker = (real_t)800; // kg/m^3
        auto ring = sys.addRigidBody(PhysicsSystem::MeshShape("res/meshes/woodpecker_ring.obj"),
                                     PhysicsSystem::RigidState(Vector3r::Zero(), Vector3r(0.0, 0.0, 0.0)),
                                     PhysicsSystem::RigidProps::withDensity(rhoRing));

        auto pecker = sys.addRigidBody(PhysicsSystem::MeshShape("res/meshes/woodpecker.obj"),
                                       PhysicsSystem::RigidState(Vector3r::Zero(), Vector3r(0.0, 0.0, -1.0)),
                                       PhysicsSystem::RigidProps::withDensity(rhoWoodpecker));

        // 3) Place hinge at ring edge (radius 5mm) in +X direction and offset woodpecker 1cm further +X
        const real_t ringEdgeRadius = (real_t)0.005; // 5mm
        const real_t woodpeckerOffset = (real_t)0.01; // 1cm from ring edge towards +X

        // 4) Add a hinge constraint between ring and woodpecker        
        const Vector3r hingeCenter = Vector3r(0, -ringEdgeRadius, 0);
        physics::JointFrame jf = physics::JointFrame(hingeCenter, Matrix33r::Identity(), std::nullopt);

        sys.addConstraint<physics::TranslationRotationConstraint>(sys.ecs(), ring, pecker, jf,
                                           /*K_trans*/ Vector3r(3000, std::numeric_limits<real_t>::infinity(), std::numeric_limits<real_t>::infinity()),
                                           /*D_trans*/ Vector3r::Zero(),
                                           /*K_rot*/   Vector3r(0.05,0.5,0.05),
                                           /*D_rot*/   Vector3r::Zero());

        // 5) Track entities for output
        sys.track(pecker, "woodpecker_body");
        sys.track(ring, "woodpecker_ring");
    }
};
