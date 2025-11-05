#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class DiscreteRodScene : public SceneBase {
public:
    DiscreteRodScene() = default;
    ~DiscreteRodScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        const size_t segments = 20;
        const real_t length = M_PI;
        const real_t segment_length = length / segments;
        const real_t width = 0.05 * length;
        const real_t height = width / 2;
        real_t mass = 12.5 / segments;

        const real_t E = 1e6;
        const real_t nu = 0.3;
        const real_t G = E / (2 * (1 + nu));
        const real_t A = width * height;
        const real_t Iy = width * std::pow(height, 3) / 12;
        const real_t Iz = std::pow(width, 3) * height / 12;
        const real_t Ip = Iy + Iz;

        Vector3r Ke(E * A, G * A, G * A);
        Vector3r Kf(G * Ip, E * Iy, E * Iz);

        PhysicsSystem::Cube shape;
        shape.halfExtents = Vector3r(segment_length / 2, width / 2, height / 2);

        Vector3r position = Vector3r::Zero();
        Quaternion4r orientation = Quaternion4r::Identity();
        Vector3r linearVelocity = Vector3r::Zero();
        Vector3r angularVelocity = Vector3r::Zero();

        entt::entity a = sys.addRigidBody(mass, position, orientation, linearVelocity, angularVelocity, shape);
        sys.makeStatic(a);

        for (index_t i=0; i < segments; ++i) {
            position += Vector3r(segment_length, 0, 0);
            // if (i == segments - 1) {
            //     mass *= 0.5;
            //     shape.halfExtents(0) *= 0.5;
            // }
            entt::entity b = sys.addRigidBody(mass, position, orientation, linearVelocity, angularVelocity, shape);
            sys.addConstraint<cardillo::physics::BeamConstraint>(sys.ecs(), a, b, Vector3r::Zero(), Vector3r::Zero(), Ke, Kf);
            a = b;
        }
    }

    // void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
    // }
};
