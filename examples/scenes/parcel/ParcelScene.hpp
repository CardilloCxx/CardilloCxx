#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class ParcelScene : public SceneBase {
public:
    const char* sceneName() const override { return "parcel"; }
    ParcelScene() = default;
    ~ParcelScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        // Create a treadmill that moves with given velocity
        real_t mass = 1e8;
        Quaternion4r orientation = Quaternion4r::Identity();
        Vector3r linearVelocity = m_linearVelocity;
        Vector3r angularVeloctiy = Vector3r::Zero();
    const real_t belt_halfLength = 5.0;
    physics::CubeShape beltShape{Vector3r(belt_halfLength, 0.5, 0.1)};
    physics::RigidState beltState; beltState.position = m_position; beltState.orientation = orientation; beltState.linearVelocity = linearVelocity; beltState.angularVelocity = angularVeloctiy;
    physics::RigidProps beltProps; beltProps.mass = mass;
    m_treatmill_entity = engine.addRigidBody(beltShape, beltState, beltProps);

        // Static walls along the conveyor (left/right)
    const real_t beltHalfX = beltShape.halfExtents.x();
    const real_t beltHalfY = beltShape.halfExtents.y();
        const real_t wallThick = (real_t)0.05;
        const real_t wallHalfZ = (real_t)0.25; // ~0.5m tall walls
        const real_t wallYoff = beltHalfY + wallThick; // place just outside belt width
        {
            physics::CubeShape wallShape{Vector3r(beltHalfX, wallThick, wallHalfZ)};
            engine.addStaticBody(wallShape, physics::RigidState(Vector3r(m_position.x(), -wallYoff, wallHalfZ)));
        }
        {
            physics::CubeShape wallShape{Vector3r(beltHalfX, wallThick, wallHalfZ)};
            engine.addStaticBody(wallShape, physics::RigidState(Vector3r(m_position.x(), +wallYoff, wallHalfZ)));
        }

        // Drop box at the end of the conveyor
        const real_t boxHalfX = (real_t)0.4;
        const real_t boxHalfY = (real_t)0.6;
        const real_t boxWallT = (real_t)0.0125;
        const real_t boxWallHalfZ = (real_t)0.4;
        const real_t boxFloorHalfZ = (real_t)0.0125;
        const real_t boxFloorZ = (real_t)-1.0; // below belt height so parcels drop in
        const real_t boxOffsetX = belt_halfLength + boxHalfX + (real_t)0.1;
        const Vector3r boxCenter(m_position.x() + boxOffsetX, 0.0, boxFloorZ);
        // Floor
        {
            physics::CubeShape floorShape{Vector3r(boxHalfX, boxHalfY, boxFloorHalfZ)};
            engine.addStaticBody(floorShape, physics::RigidState(Vector3r(boxCenter.x(), boxCenter.y(), boxFloorZ)));
        }
        // Back wall
        {
            physics::CubeShape backShape{Vector3r(boxWallT, boxHalfY, 2 * boxWallHalfZ)};
            engine.addStaticBody(backShape, physics::RigidState(Vector3r(boxCenter.x() + boxHalfX + boxWallT, boxCenter.y(), boxCenter.z() + 2 * boxWallHalfZ)));
        }
        // Front wall
        {
            physics::CubeShape frontShape{Vector3r(boxWallT, boxHalfY, boxWallHalfZ)};
            engine.addStaticBody(frontShape, physics::RigidState(Vector3r(boxCenter.x() - (boxHalfX + boxWallT), boxCenter.y(), boxCenter.z() + boxWallHalfZ)));
        }
        // Side walls
        {
            physics::CubeShape sideLShape{Vector3r(boxHalfX + boxWallT, boxWallT, boxWallHalfZ)};
            engine.addStaticBody(sideLShape, physics::RigidState(Vector3r(boxCenter.x(), boxCenter.y() - (boxHalfY + boxWallT), boxCenter.z() + boxWallHalfZ)));
        }
        {
            physics::CubeShape sideRShape{Vector3r(boxHalfX + boxWallT, boxWallT, boxWallHalfZ)};
            engine.addStaticBody(sideRShape, physics::RigidState(Vector3r(boxCenter.x(), boxCenter.y() + (boxHalfY + boxWallT), boxCenter.z() + boxWallHalfZ)));
        }
    }

    void maybeSpawnParcel(physics::PhysicsEngine& engine, real_t t) {
        using namespace cardillo;

        if (t > m_t0) {
            m_t0 += 0.5;
            real_t mass = 1.0;
            // Vector3r position = Vector3r::Random();
            Vector3r position(0.0, 0.0, 1.0);
            Quaternion4r orientation = Quaternion4r::UnitRandom();
            Vector3r linearVelocity = Vector3r::Zero();
            Vector3r angularVelocity = Vector3r::Zero();
            physics::CubeShape parcelShape{0.1 * (Vector3r::Random() + Vector3r::Constant(1.2)) / 2.2};
            physics::RigidState st; st.position = position; st.orientation = orientation; st.linearVelocity = linearVelocity; st.angularVelocity = angularVelocity; physics::RigidProps pr; pr.mass = mass; engine.addRigidBody(parcelShape, st, pr);
        }
    }

    void resetPosition(physics::PhysicsEngine& engine) {
        engine.setPosition(m_treatmill_entity, m_position);
        engine.setOrientation(m_treatmill_entity, Quaternion4r::Identity());
        engine.setLinearVelocity(m_treatmill_entity, m_linearVelocity);
        engine.setAngularVelocity(m_treatmill_entity, Vector3r::Zero());
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        maybeSpawnParcel(engine, t);
        if (m_treatmill_entity != entt::null) {
            resetPosition(engine);
        }
    }
private:
    entt::entity m_treatmill_entity = entt::null;
    real_t m_t0 = 0.0;
    Vector3r m_position = Vector3r(3.5, 0.0, -0.1);
    real_t m_v0 = 4.0;
    Vector3r m_linearVelocity = Vector3r(m_v0, 0.0, 0.0);
};
