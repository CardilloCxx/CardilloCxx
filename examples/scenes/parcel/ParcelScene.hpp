#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Parcel scene.
class ParcelScene : public SceneBase {
public:
    ParcelScene() = default;
    ~ParcelScene() = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Create a treadmill that moves with given velocity
        real_t mass = 1e8;
        Quaternion4r orientation = Quaternion4r::Identity();
        Vector3r linearVelocity = m_linearVelocity;
        Vector3r angularVeloctiy = Vector3r::Zero();
        PhysicsSystem::Cube shape;
        shape.center = m_position;
        const real_t belt_halfLength = 5.0;
        shape.halfExtents = Vector3r(belt_halfLength, 0.5, 0.1);
        m_treatmill_entity = sys.addRigidBody(mass, m_position, orientation, linearVelocity, angularVeloctiy, shape);

        // Static walls along the conveyor (left/right)
        const real_t beltHalfX = shape.halfExtents.x();
        const real_t beltHalfY = shape.halfExtents.y();
        const real_t wallThick = (real_t)0.05;
        const real_t wallHalfZ = (real_t)0.25; // ~0.5m tall walls
        const real_t wallYoff = beltHalfY + wallThick; // place just outside belt width
        {
            PhysicsSystem::Cube wallL;
            wallL.center = Vector3r(m_position.x(), -wallYoff, wallHalfZ);
            wallL.halfExtents = Vector3r(beltHalfX, wallThick, wallHalfZ);
            wallL.q = Quaternion4r::Identity();
            sys.addObstacleBody(wallL);
        }
        {
            PhysicsSystem::Cube wallR;
            wallR.center = Vector3r(m_position.x(), +wallYoff, wallHalfZ);
            wallR.halfExtents = Vector3r(beltHalfX, wallThick, wallHalfZ);
            wallR.q = Quaternion4r::Identity();
            sys.addObstacleBody(wallR);
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
            PhysicsSystem::Cube floor;
            floor.center = Vector3r(boxCenter.x(), boxCenter.y(), boxFloorZ);
            floor.halfExtents = Vector3r(boxHalfX, boxHalfY, boxFloorHalfZ);
            floor.q = Quaternion4r::Identity();
            sys.addObstacleBody(floor);
        }
        // Back wall
        {
            PhysicsSystem::Cube back;
            back.center = Vector3r(boxCenter.x() + boxHalfX + boxWallT, boxCenter.y(), boxCenter.z() + 2 * boxWallHalfZ);
            back.halfExtents = Vector3r(boxWallT, boxHalfY, 2 * boxWallHalfZ);
            back.q = Quaternion4r::Identity();
            sys.addObstacleBody(back);
        }
        // Front wall
        {
            PhysicsSystem::Cube front;
            front.center = Vector3r(boxCenter.x() - (boxHalfX + boxWallT), boxCenter.y(), boxCenter.z() + boxWallHalfZ);
            front.halfExtents = Vector3r(boxWallT, boxHalfY, boxWallHalfZ);
            front.q = Quaternion4r::Identity();
            sys.addObstacleBody(front);
        }
        // Side walls
        {
            PhysicsSystem::Cube sideL;
            sideL.center = Vector3r(boxCenter.x(), boxCenter.y() - (boxHalfY + boxWallT), boxCenter.z() + boxWallHalfZ);
            sideL.halfExtents = Vector3r(boxHalfX + boxWallT, boxWallT, boxWallHalfZ);
            sideL.q = Quaternion4r::Identity();
            sys.addObstacleBody(sideL);
        }
        {
            PhysicsSystem::Cube sideR;
            sideR.center = Vector3r(boxCenter.x(), boxCenter.y() + (boxHalfY + boxWallT), boxCenter.z() + boxWallHalfZ);
            sideR.halfExtents = Vector3r(boxHalfX + boxWallT, boxWallT, boxWallHalfZ);
            sideR.q = Quaternion4r::Identity();
            sys.addObstacleBody(sideR);
        }
    }

    void maybeSpawnParcel(cardillo::PhysicsSystem& sys, real_t t) {
        using namespace cardillo;

        if (t > m_t0) {
            m_t0 += 0.5;
            real_t mass = 1.0;
            // Vector3r position = Vector3r::Random();
            Vector3r position(0.0, 0.0, 1.0);
            Quaternion4r orientation = Quaternion4r::UnitRandom();
            Vector3r linearVelocity = Vector3r::Zero();
            Vector3r angularVelocity = Vector3r::Zero();
            PhysicsSystem::Cube shape;
            shape.center = position;
            shape.halfExtents = 0.1 * (Vector3r::Random() + Vector3r::Constant(1.2)) / 2.2;
            shape.q = orientation;
            sys.addRigidBody(mass, position, orientation, linearVelocity, angularVelocity, shape);
        }
    }

    void resetPosition(cardillo::PhysicsSystem& sys) {
        sys.setPosition(m_treatmill_entity, m_position);
        sys.setOrientation(m_treatmill_entity, Quaternion4r::Identity());
        sys.setLinearVelocity(m_treatmill_entity, m_linearVelocity);
        sys.setAngularVelocity(m_treatmill_entity, Vector3r::Zero());
    }

    void updateScene(cardillo::PhysicsSystem& sys, real_t t, real_t /*dt*/) override {
        maybeSpawnParcel(sys, t);
        if (m_treatmill_entity != entt::null) {
            resetPosition(sys);
        }
    }
private:
    entt::entity m_treatmill_entity = entt::null;
    real_t m_t0 = 0.0;
    Vector3r m_position = Vector3r(3.5, 0.0, -0.1);
    real_t m_v0 = 4.0;
    Vector3r m_linearVelocity = Vector3r(m_v0, 0.0, 0.0);
};
