#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

using namespace cardillo;

class BeamTestScene : public SceneBase {
public:
    BeamTestScene() = default;
    ~BeamTestScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        size_t segments = 50;
        real_t beamLength = 2.0;
        real_t segmentLength = beamLength / real_t(segments);
        real_t beamWidth = 0.1;
        real_t beamDensity = 500.0; // kg/m^3

        PhysicsSystem::Cube cshape; cshape.halfExtents = Vector3r(segmentLength/2, beamWidth/2, beamWidth/2);
        Quaternion4r q = Quaternion4r::Identity();
        Vector3r v0 = Vector3r::Zero();
        Vector3r w0 = Vector3r::Zero();
        real_t segmentVolume = segmentLength * beamWidth * beamWidth;
        real_t segmentMass = beamDensity * segmentVolume;
        entt::entity A = sys.addRigidBody(segmentMass, Vector3r(0.0 - beamLength/2, 0, 0.0), q, v0, w0, cshape);
        sys.makeStatic(A);

        // Attachment at cube centers
        Vector3r Ke(1e4, 1e4, 1e4); // stretch/shear stiffnesses
        Vector3r Kf(5e3, 5e3, 5e3); // torsion/bend stiffnesses
        Vector3r De(0, 0, 0); // stretch/shear damping
        Vector3r Df(0, 0, 0); // torsion/bend damping

        for (size_t i = 1; i < segments; ++i) {
            real_t alpha = real_t(i) / real_t(segments);
            Vector3r pos = Vector3r(-beamLength/2 + alpha * beamLength, 0.0, 0.0);
            entt::entity B = sys.addRigidBody(segmentMass, pos, q, v0, w0, cshape);
            sys.addConstraint<physics::BeamConstraint>(sys.ecs(), A, B, Vector3r::Zero(), Vector3r::Zero(), Ke, Kf, De, Df);
            A = B; // chain
            m_beamEnd = B;
        }
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t /*t*/, real_t /*dt*/) override {
        sys.applyForce(m_beamEnd, Vector3r(0, 0, 0), Vector3r(0, 0, 100));
    }

private:
    entt::entity m_beamEnd{entt::null};
};
