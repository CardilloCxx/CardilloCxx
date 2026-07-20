#pragma once

#include "../SceneBase.hpp"
#include "misc/spline.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace cardillo;

class ChainFountainScene : public SceneBase {
public:
    const char* sceneName() const override { return "chain_fountain"; }
    ChainFountainScene() = default;
    ~ChainFountainScene() = default;

    Vector3r m_containerDims = Vector3r(0.04, 0.04, 0.025); 
    real_t m_capRadius = 0.0015;
    real_t m_capLength = 0.0000001;
    real_t m_distanceBetweenCaps = 0.0002;
    real_t m_spacing = 2.0 * m_capRadius + m_distanceBetweenCaps;        
    int m_numLayers = 10;         
    real_t m_layerSpacing = m_capRadius * 3.0;
    real_t m_dropHeight = 20.0;
    real_t m_settlingTime = 0.25;

    real_t m_yeetUpTime = 0.2;
    real_t m_yeetArcTime = 0.1;
    real_t m_yeetDownTime = 0.2;
    real_t m_yeetForceMag = 0.1;
    Vector3r m_yeetHorizontalDir = Vector3r(1.0, 0.0, 0.0);
    real_t m_yeetStartAngleDeg = 90.0;
    real_t m_yeetEndAngleDeg = -90.0;

    std::vector<entt::entity> m_chainEntities;
    bool m_yeeted = false;
    real_t m_yeetStartTime = -1.0;

private:
    real_t findNextAlpha(const misc::CatmullRomSpline& spline, real_t startAlpha, const Vector3r& startP, real_t targetDist) {
        real_t a = startAlpha;
        const real_t step = 0.0005;
        
        while (a < 1.0) {
            a += step;
            bool clamped = false;
            
            if (a >= 1.0) {
                a = 1.0;
                clamped = true;
            }
            
            Vector3r p = spline.sample(a).position;
            real_t d = (p - startP).norm();
            
            if (d >= targetDist) {
                real_t safePrevAlpha = std::max(static_cast<real_t>(0.0), a - step);
                Vector3r prevP = spline.sample(safePrevAlpha).position;
                real_t d1 = (prevP - startP).norm();
                real_t f = (d > d1) ? (targetDist - d1) / (d - d1) : 0.0;
                return (a - step) + f * step;
            }
            
            if (clamped) {
                break;
            }
        }
        return 1.1;
    }

public:
    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        engine.setGravity(Vector3r(0, 0, -9.81));

        physics::CubeShape groundShape(Vector3r(10.0, 10.0, 0.1));
        engine.addStaticBody(groundShape, physics::RigidState(Vector3r(0.0, 0.0, -0.1)));

        const real_t wallThickness = 0.005;
        const real_t halfThick = wallThickness * 0.5;
        const real_t cx = m_containerDims.x() * 0.5;
        const real_t cy = m_containerDims.y() * 0.5;
        const real_t cz = m_containerDims.z() * 0.5;
        const Vector3r center(0.0, 0.0, m_dropHeight + cz);

        engine.addStaticBody(physics::CubeShape(Vector3r(cx + wallThickness, cy + wallThickness, halfThick)),
            physics::RigidState(center + Vector3r(0.0, 0.0, -cz - halfThick)));
        engine.addStaticBody(physics::CubeShape(Vector3r(halfThick, cy + wallThickness, cz)),
            physics::RigidState(center + Vector3r(-cx - halfThick, 0.0, 0.0)));
        engine.addStaticBody(physics::CubeShape(Vector3r(halfThick, cy + wallThickness, cz)),
            physics::RigidState(center + Vector3r(cx + halfThick, 0.0, 0.0)));
        engine.addStaticBody(physics::CubeShape(Vector3r(cx, halfThick, cz)),
            physics::RigidState(center + Vector3r(0.0, -cy - halfThick, 0.0)));
        engine.addStaticBody(physics::CubeShape(Vector3r(cx, halfThick, cz)),
            physics::RigidState(center + Vector3r(0.0, cy + halfThick, 0.0)));

        std::vector<Vector3r> splinePoints;
        const real_t padding = m_capRadius;
        const real_t usableX = cx - padding;
        const real_t usableY = cy - padding;
        const real_t rowSpacing = m_capRadius * 3.0;

        int numRowsPerLayer = std::max(1, static_cast<int>(std::floor(2.0 * usableY / rowSpacing)) + 1);
        real_t actualRowSpacing = (2.0 * usableY) / std::max(1, numRowsPerLayer - 1);

        real_t currZ = m_dropHeight + padding;
        real_t currY = -usableY;
        real_t currX = -usableX;
        
        bool sweepXPositive = true;
        bool sweepYPositive = true;

        for (int layer = 0; layer < m_numLayers; ++layer) {
            for (int row = 0; row < numRowsPerLayer; ++row) {
                
                real_t startX = sweepXPositive ? -usableX : usableX;
                real_t endX = sweepXPositive ? usableX : -usableX;
                for (int s = 0; s <= 10; ++s) {
                    real_t t = static_cast<real_t>(s) / 10.0;
                    splinePoints.push_back(Vector3r(startX + t * (endX - startX), currY, currZ));
                }
                currX = endX;
                sweepXPositive = !sweepXPositive;
                
                if (row < numRowsPerLayer - 1) {
                    real_t nextY = currY + (sweepYPositive ? actualRowSpacing : -actualRowSpacing);
                    for (int s = 1; s <= 5; ++s) {
                        real_t t = static_cast<real_t>(s) / 5.0;
                        splinePoints.push_back(Vector3r(currX, currY + t * (nextY - currY), currZ));
                    }
                    currY = nextY;
                }
            }
            
            if (layer < m_numLayers - 1) {
                real_t nextZ = currZ + m_layerSpacing;
                for (int s = 1; s <= 5; ++s) {
                    real_t t = static_cast<real_t>(s) / 5.0;
                    splinePoints.push_back(Vector3r(currX, currY, currZ + t * (nextZ - currZ)));
                }
                currZ = nextZ;
                sweepYPositive = !sweepYPositive;
            }
        }

        CatmullRomSpline chainSpline(splinePoints, false);

        physics::SphereShape sphereShape(m_capRadius);
        
        physics::RigidProps capProps = physics::RigidProps::withDensity(8000.0);
        capProps.restitution_normal = 0.1;
        capProps.restitution_tangential = 0.05;
        real_t currentAlpha = 0.0;
        
        while (currentAlpha <= 1.0) {
            Vector3r cap1Pos = chainSpline.sample(currentAlpha).position;
            
            real_t nextAlpha = findNextAlpha(chainSpline, currentAlpha, cap1Pos, m_capLength);
            if (nextAlpha > 1.0) break;
            Vector3r cap2Pos = chainSpline.sample(nextAlpha).position;

            Vector3r center = (cap1Pos + cap2Pos) * 0.5;
            
            Vector3r dir = (cap2Pos - cap1Pos).normalized();
            Quaternion4r rot = Eigen::Quaternion<real_t>::FromTwoVectors(Vector3r::UnitZ(), dir);
            
            auto capEntity = engine.addRigidBody(sphereShape, physics::RigidState(center, rot), capProps);
            m_chainEntities.push_back(capEntity);
            
            currentAlpha = findNextAlpha(chainSpline, nextAlpha, cap2Pos, m_spacing);
        }
        engine.track(m_chainEntities.back(), "chain_end");
        size_t numCaps = m_chainEntities.size();
        real_t chainLength = (numCaps - 1) * (m_spacing + m_capLength);
        std::cout << "ChainFountainScene: Created " << numCaps << " capsules, total chain length = " << chainLength << std::endl;


        const Vector3r localTop(0.0, 0.0, m_capLength * 0.5 + m_capRadius);
        const Vector3r localBottom(0.0, 0.0, -m_capLength * 0.5 - m_capRadius);

        for (size_t i = 1; i < m_chainEntities.size(); ++i) {
            auto capPrev = m_chainEntities[i - 1];
            auto capCurr = m_chainEntities[i];

            engine.addLinearDistanceConstraint(capPrev, capCurr, localTop, localBottom, std::numeric_limits<real_t>::infinity(), 0.0, m_distanceBetweenCaps);
        }
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {

        if (t >= m_settlingTime && !m_yeeted && !m_chainEntities.empty()) {
            m_yeeted = true;
            m_yeetStartTime = t;

//             auto& registry = engine.world().ecs();
//             for (auto& cap : m_chainEntities) {
//                 if (registry.any_of<cardillo::physics::RigidProps>(cap)) {
//                     auto& capProps = registry.get<cardillo::physics::RigidProps>(cap);
// 
//                     capProps.restitution_normal = 0.99;
//                     capProps.restitution_tangential = 0.1;
//                 }
//             }
        }

        if (m_yeeted && !m_chainEntities.empty()) {
            const real_t elapsed = t - m_yeetStartTime;

            if (elapsed >= 0.0 && elapsed <= m_yeetUpTime + m_yeetArcTime + m_yeetDownTime) {
                const real_t s = elapsed / (m_yeetUpTime + m_yeetArcTime + m_yeetDownTime);
                
                const real_t upS = std::clamp(elapsed / m_yeetUpTime, 0.0, 1.0);
                const real_t arcS = std::clamp((elapsed - m_yeetUpTime) / m_yeetArcTime, 0.0, 1.0);
                const real_t downS = std::clamp((elapsed - m_yeetUpTime - m_yeetArcTime) / m_yeetDownTime, 0.0, 1.0);
                const real_t envelope = upS + downS;

                constexpr real_t kPi = static_cast<real_t>(3.14159265358979323846);
                const real_t angleDeg = m_yeetStartAngleDeg + (m_yeetEndAngleDeg - m_yeetStartAngleDeg) * arcS;
                const real_t angleRad = angleDeg * kPi / static_cast<real_t>(180.0);

                Vector3r dir = m_yeetHorizontalDir.normalized() * std::cos(angleRad) + Vector3r::UnitZ() * std::sin(angleRad);

                auto endCapsule = m_chainEntities.back();
                engine.applyForce(endCapsule, dir * (m_yeetForceMag * envelope), Vector3r::Zero());
            }
        }
    }
};