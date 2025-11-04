#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class ParcelScene {
public:
    ParcelScene() = default;
    ~ParcelScene() = default;

    void populate(cardillo::PhysicsSystem& sys) {
        using namespace cardillo;

        // Create a treadmill that moves with given velocity
        real_t mass = 1e8;
        Quaternion4r orientation = Quaternion4r::Identity();
        Vector3r linearVelocity = m_linearVelocity;
        Vector3r angularVeloctiy = Vector3r::Zero();
        PhysicsSystem::Cube shape;
        shape.center = m_position;
        shape.halfExtents = Vector3r(5.0, 0.5, 0.1);
        m_treatmill_entity = sys.addRigidBody(mass, m_position, orientation, linearVelocity, angularVeloctiy, shape);
    }

    void maybeSpawnParcel(cardillo::PhysicsSystem& sys, real_t t) {
        using namespace cardillo;

        if (t > m_t0) {
        // if (t > 1.0 && t < 1.002) {
            m_t0 += 0.5;
            real_t mass = 1.0;
            // // Vector3r position = Vector3r::Random();
            // Vector3r position = Vector3r::Zero();
            // position(0) = 0.0;
            // position(1) = 0.0;
            // // position(2) += 1.0;
            // position(2) = 5.0 + t;
            Vector3r position(0.0, 0.0, 1.0);
            Quaternion4r orientation = Quaternion4r::UnitRandom();
            Vector3r linearVelocity = Vector3r::Zero();
            Vector3r angularVeloctiy = Vector3r::Zero();
            PhysicsSystem::Cube shape;
            shape.center = position;
            shape.halfExtents = 0.1 * (Vector3r::Random() + Vector3r::Constant(1.0)) / 2.0;
            shape.q = orientation;
            sys.addRigidBody(mass, position, orientation, linearVelocity, angularVeloctiy, shape);
        }
    }

    void resetPosition(cardillo::PhysicsSystem& sys) {
        sys.ecs().get<PhysicsSystem::C_Position3>(m_treatmill_entity).value = m_position;
        sys.ecs().get<PhysicsSystem::C_Orientation>(m_treatmill_entity).q = Quaternion4r::Identity();
        sys.ecs().get<PhysicsSystem::C_LinearVelocity3>(m_treatmill_entity).value = m_linearVelocity;
    }
private:
    entt::entity m_treatmill_entity = entt::null;
    real_t m_t0 = 0.0;
    Vector3r m_position = Vector3r(3.5, 0.0, -0.1);
    real_t m_v0 = 1.0;
    Vector3r m_linearVelocity = Vector3r(m_v0, 0.0, 0.0);
};
