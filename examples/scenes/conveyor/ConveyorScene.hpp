#pragma once

#include "../SceneBase.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace cardillo;

// https://www.youtube.com/watch?v=6PzVqoEKuoM

class ConveyorScene : public SceneBase {
   public:
    const char* sceneName() const override { return "conveyor"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        m_engine = &engine;
        engine.setGravity(Vector3r(0, 0, -9.81));

        auto makeState = [](const Vector3r& p, real_t yaw, real_t pitch) {
            physics::RigidState s;
            s.position = p;
            const Quaternion4r qYaw(Eigen::AngleAxis<real_t>(yaw, Vector3r::UnitZ()));
            const Quaternion4r qPitch(Eigen::AngleAxis<real_t>(pitch, Vector3r::UnitY()));
            s.orientation = qYaw * qPitch;
            return s;
        };

        const real_t speed = (real_t)0.8;
        const real_t w = (real_t)0.20;
        const real_t h = (real_t)0.03;
        const size_t n = (size_t)72;

        // Top flat orthogonal conveyor (runs along +Y).
        const real_t topLen = (real_t)0.70;
        const real_t topRadius = (real_t)0.5 * h;
        const real_t topRunZ = (real_t)0.20;
        physics::RigidState topState = makeState(Vector3r((real_t)0.0, (real_t)0.0, topRunZ - topRadius), (real_t)M_PI_2, (real_t)0.0);
        createConveyor(topState, topLen, w, h, n, speed, +1);

        // Two symmetric ramps feeding the top conveyor at the back side edges.
        const real_t rampLen = (real_t)0.55;
        const real_t rampRadius = (real_t)0.5 * h;
        const real_t rampAngle = (real_t)(20.0 * M_PI / 180.0);
        const real_t yBackTop = (real_t)-0.5 * topLen;
        const real_t yRampCenter = yBackTop + (real_t)0.5 * w;  // back edge of side ramps flush with top-back plane

        const Vector3r leftJoin((real_t)-0.5 * w, yRampCenter, topRunZ);
        const Vector3r rightJoin((real_t) + 0.5 * w, yRampCenter, topRunZ);

        const real_t leftYaw = (real_t)0.0;
        const real_t rightYaw = (real_t)0.0;

        physics::RigidState leftRampState = makeState(Vector3r((real_t)-0.30, yRampCenter, topRunZ - (real_t)0.11), leftYaw, -rampAngle);
        const Matrix33r RL = leftRampState.orientation.toRotationMatrix();
        leftRampState.position =
            leftJoin - RL * Vector3r((real_t)0.5 * rampLen, (real_t)0.0, rampRadius) + Vector3r(-0.015, 0.0, 0.0);  // small extra forward offset to prevent initial interpenetration with top conveyor
        createConveyor(leftRampState, rampLen, w, h, n, speed, +1);

        physics::RigidState rightRampState = makeState(Vector3r((real_t)0.30, yRampCenter, topRunZ - (real_t)0.11), rightYaw, +rampAngle);
        const Matrix33r RR = rightRampState.orientation.toRotationMatrix();
        rightRampState.position = rightJoin - RR * Vector3r((real_t)-0.5 * rampLen, (real_t)0.0, rampRadius) +
                                  Vector3r(+0.015, 0.0, 0.0);  // small extra forward offset to prevent initial interpenetration with top conveyor
        createConveyor(rightRampState, rampLen, w, h, n, speed, -1);

        // Two additional flat feeders, one feeding each diagonal ramp at its lower edge.
        const real_t feedLen = (real_t)0.42;

        // Formal world-geometry anchors for feeder attachment at ramp lower bottom endpoints.
        // Left ramp join used local (+L/2, z=+r), so lower opposite bottom endpoint is local (-L/2, z=0).
        const Vector3r leftRampBottom = leftRampState.position + RL * Vector3r((real_t)-0.5 * rampLen, (real_t)0.0, (real_t)0.0) +
                                        Vector3r(-0.025, 0.0, 0.0);  // small extra rearward offset to prevent initial interpenetration with ramp
        // Right ramp join used local (-L/2, z=+r), so lower opposite bottom endpoint is local (+L/2, z=0).
        const Vector3r rightRampBottom = rightRampState.position + RR * Vector3r((real_t) + 0.5 * rampLen, (real_t)0.0, (real_t)0.0) +
                                         Vector3r(+0.025, 0.0, 0.0);  // small extra rearward offset to prevent initial interpenetration with ramp

        physics::RigidState leftFeedState = makeState(Vector3r::Zero(), (real_t)0.0, (real_t)0.0);
        leftFeedState.position = leftRampBottom - Vector3r((real_t)0.5 * feedLen, (real_t)0.0, (real_t)0.0);
        createConveyor(leftFeedState, feedLen, w, h, n, speed * (real_t)0.9, +1);

        physics::RigidState rightFeedState = makeState(Vector3r::Zero(), (real_t)M_PI, (real_t)0.0);
        rightFeedState.position = rightRampBottom + Vector3r((real_t)0.5 * feedLen, (real_t)0.0, (real_t)0.0);
        createConveyor(rightFeedState, feedLen, w, h, n, speed * (real_t)0.9, +1);

        // Spawn positions centered above each feeder belt.
        const real_t spawnLift = (real_t)0.08;
        m_leftFeederSpawn = leftFeedState.position + Vector3r((real_t)0.0, (real_t)0.0, (real_t)0.5 * h + spawnLift);
        m_rightFeederSpawn = rightFeedState.position + Vector3r((real_t)0.0, (real_t)0.0, (real_t)0.5 * h + spawnLift);

        // Orthogonal take-away conveyor at the end of top belt (runs along +X).
        const real_t outLen = w;
        const real_t outRadius = (real_t)0.5 * h;
        const real_t yFrontTop = (real_t)0.5 * topLen;
        const real_t yOutCenter = yFrontTop + (real_t)0.5 * w + 0.005;  // back edge of out belt flush with top-front plane
        physics::RigidState outState = makeState(Vector3r((real_t)0.0, yOutCenter, topRunZ - outRadius - 0.005), (real_t)0.0, (real_t)0.0);
        createConveyor(outState, outLen, w, h, n / 2, speed / 2.0, +1);

        // Slightly downward roller array at the outfeed end.
        const real_t rollerRadius = (real_t)0.012;
        const real_t rollerSpacing = (real_t)2.2 * rollerRadius;
        const real_t rollerRunLen = (real_t)4.0 * outLen;
        const Matrix33r ROut = outState.orientation.toRotationMatrix();
        const Vector3r outX = ROut.col(0);
        const Vector3r outZ = ROut.col(2);
        const Vector3r outTopFront = outState.position + outX * ((real_t)0.5 * outLen) + outZ * ((real_t)0.5 * h);
        const real_t rollerHalfLength = (real_t)0.49 * w;

        // Build rollers from explicit start/end coordinates.
        const Vector3r rollerStart = outTopFront + 3.0 * outX * rollerRadius - outZ * rollerRadius;
        const real_t rollerSlope = std::tan((real_t)(12.0 * M_PI / 180.0));
        const Vector3r rollerEnd = rollerStart + outX * rollerRunLen - outZ * (rollerSlope * rollerRunLen);
        createRollers(rollerStart, rollerEnd, rollerSpacing, rollerRadius, rollerHalfLength, (real_t)0.000005);
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        // Spawn new packages at regular intervals with random shape and initial velocity.
        const real_t spawnInterval = (real_t)0.5;
        if (t >= m_nextSpawnTime) {
            m_nextSpawnTime += spawnInterval;
            const bool spawnLeft = (std::rand() % 2) == 0;
            const Vector3r spawnPos = spawnLeft ? m_leftFeederSpawn : m_rightFeederSpawn;
            const real_t sizeX = (real_t)0.02 + (real_t)0.08 * ((real_t)std::rand() / (real_t)RAND_MAX);
            const real_t sizeY = (real_t)0.02 + (real_t)0.08 * ((real_t)std::rand() / (real_t)RAND_MAX);
            const real_t sizeZ = (real_t)0.02 + (real_t)0.08 * ((real_t)std::rand() / (real_t)RAND_MAX);
            const Vector3r halfExtents(sizeX * (real_t)0.5, sizeY * (real_t)0.5, sizeZ * (real_t)0.5);
            physics::RigidState state(spawnPos);
            state.linearVelocity = Vector3r::Random() * (real_t)0.2;
            state.angularVelocity = Vector3r::Random() * (real_t)4.0;
            engine.addRigidBody(physics::CubeShape(halfExtents), state, physics::RigidProps::withDensity((real_t)500.0));
        }
    }

   private:
    std::vector<entt::entity> createRollers(const Vector3r& start, const Vector3r& end, real_t spacing, real_t radius, real_t halfLength, real_t hingeDamping) {
        std::vector<entt::entity> rollers;
        if (!m_engine || spacing <= (real_t)0.0) return rollers;

        const Vector3r line = end - start;
        const real_t lineLen = line.norm();
        if (lineLen <= (real_t)1e-9) return rollers;

        const Vector3r lineDir = line / lineLen;
        const size_t count = (size_t)std::floor(lineLen / spacing) + 1;

        const Vector3r axis = Vector3r::UnitY();
        const Quaternion4r qCylAxis(Eigen::AngleAxis<real_t>(-(real_t)M_PI_2, Vector3r::UnitX()));

        for (size_t i = 0; i < count; ++i) {
            const Vector3r pWorld = start + lineDir * ((real_t)i * spacing);

            physics::RigidState rs;
            rs.position = pWorld;
            rs.orientation = qCylAxis;

            physics::RigidProps props((real_t)0.05);
            props.friction = (real_t)1.2;
            entt::entity roller = m_engine->addRigidBody(physics::CylinderShape(radius, halfLength), rs, props);
            rollers.push_back(roller);

            physics::JointFrame jf = physics::JointFrame::fromAxis(pWorld, axis);
            m_engine->addHingeConstraint(roller, entt::null, jf, (real_t)0.0, hingeDamping);
        }

        return rollers;
    }

   private:
    std::vector<entt::entity> createConveyor(const physics::RigidState& state, real_t length, real_t width, real_t height, size_t segments, real_t speed = (real_t)0.55, int direction = +1) {
        std::vector<entt::entity> elements;
        if (!m_engine || segments < 4) return elements;

        const real_t radius = std::max((real_t)0.002, (real_t)0.5 * height);
        const real_t straight = std::max((real_t)0.05, length);
        const real_t loopLength = (real_t)2.0 * straight + (real_t)2.0 * M_PI * radius;
        const real_t ds = loopLength / (real_t)segments;
        const real_t signedSpeed = ((direction >= 0) ? (real_t)1 : (real_t)-1) * std::abs(speed);

        const real_t slatPitch = ds;
        const real_t slatThickness = std::max((real_t)0.001, (real_t)0.18 * radius);
        const Vector3r slatHalfExtents((real_t)0.5 * slatPitch, (real_t)0.5 * width, (real_t)0.5 * slatThickness);

        const Matrix33r R0 = state.orientation.toRotationMatrix();
        const Vector3r p0 = state.position;

        auto wrapS = [=](real_t s) {
            real_t out = std::fmod(s, loopLength);
            if (out < (real_t)0.0) out += loopLength;
            return out;
        };

        auto sampleCenterline = [=](real_t s) {
            real_t x = 0;
            real_t z = 0;

            if (s < straight) {
                const real_t u = s / straight;
                x = (-(real_t)0.5 + u) * straight;
                z = radius;
            } else if (s < straight + M_PI * radius) {
                const real_t u = (s - straight) / (M_PI * radius);
                const real_t theta = (real_t)M_PI_2 - u * (real_t)M_PI;
                x = (real_t)0.5 * straight + radius * std::cos(theta);
                z = radius * std::sin(theta);
            } else if (s < (real_t)2.0 * straight + M_PI * radius) {
                const real_t u = (s - straight - M_PI * radius) / straight;
                x = ((real_t)0.5 - u) * straight;
                z = -radius;
            } else {
                const real_t u = (s - ((real_t)2.0 * straight + M_PI * radius)) / (M_PI * radius);
                const real_t theta = -(real_t)M_PI_2 - u * (real_t)M_PI;
                x = -(real_t)0.5 * straight + radius * std::cos(theta);
                z = radius * std::sin(theta);
            }

            return Vector3r(x, (real_t)0.0, z);
        };

        auto sampleSlatPose = [=](real_t sMid) {
            const real_t halfLen = slatHalfExtents.x();
            const Vector3r pBack = sampleCenterline(wrapS(sMid - halfLen));
            const Vector3r pFront = sampleCenterline(wrapS(sMid + halfLen));

            Vector3r xAxis = pFront - pBack;
            if (xAxis.squaredNorm() < (real_t)1e-16) {
                const Vector3r pA = sampleCenterline(wrapS(sMid - (real_t)1e-4));
                const Vector3r pB = sampleCenterline(wrapS(sMid + (real_t)1e-4));
                xAxis = pB - pA;
            }
            xAxis.normalize();

            Vector3r yAxis = Vector3r::UnitY();
            Vector3r zAxis = xAxis.cross(yAxis).normalized();
            yAxis = zAxis.cross(xAxis).normalized();

            Matrix33r R = Matrix33r::Identity();
            R.col(0) = xAxis;
            R.col(1) = yAxis;
            R.col(2) = zAxis;

            const Vector3r pCenter = ((real_t)0.5) * (pBack + pFront);
            return std::pair<Vector3r, Quaternion4r>(pCenter, Quaternion4r(R));
        };

        physics::RigidProps slatProps((real_t)0.3);
        slatProps.friction = (real_t)1.0;

        for (size_t i = 0; i < segments; ++i) {
            const real_t phase = (real_t)i / (real_t)segments;

            const real_t s0 = phase * loopLength;
            const auto [pLocal0, qLocal0] = sampleSlatPose(s0);

            physics::RigidState rs;
            rs.position = p0 + R0 * pLocal0;
            rs.orientation = state.orientation * qLocal0;

            const entt::entity e = m_engine->addRigidBody(physics::CubeShape(slatHalfExtents), rs, slatProps);
            elements.push_back(e);

            m_engine->addTrajectory(e, std::make_optional<std::function<TrajectoryPose(real_t)>>([=](real_t t) -> TrajectoryPose {
                                        const real_t s = wrapS(s0 + signedSpeed * t);
                                        const auto [pLocal, qLocal] = sampleSlatPose(s);
                                        return TrajectoryPose{p0 + R0 * pLocal, state.orientation * qLocal};
                                    }),
                                    std::nullopt);
        }

        return elements;
    }

   private:
    cardillo::physics::PhysicsEngine* m_engine{nullptr};
    real_t m_nextSpawnTime{0};
    Vector3r m_leftFeederSpawn{Vector3r::Zero()};
    Vector3r m_rightFeederSpawn{Vector3r::Zero()};
};
