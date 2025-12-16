#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace cardillo;

// Cardhouse scene: stack real-world playing cards into layered A-frames with
// flat roof cards. The geometry is deterministic and uses realistic card
// dimensions so the stack can be reproduced across runs.
class CardhouseScene : public SceneBase {
public:
    const char* sceneName() const override { return "cardhouse"; }
    CardhouseScene() = default;
    ~CardhouseScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Ground plane
        PhysicsSystem::CubeShape groundShape{Vector3r(2.0, 2.0, 0.1)};
        PhysicsSystem::RigidState groundState; groundState.position = Vector3r(0.0, 0.0, -0.1);
        sys.addStaticBody(groundShape, groundState);

        // Real card dimensions (meters): long=88.9mm along x, short=63.5mm along y, thickness ~0.3 mm along z
        const Vector3r halfExtents((real_t)0.0889 * (real_t)0.5, (real_t)0.0635 * (real_t)0.5, (real_t)0.0003 * (real_t)0.5);
        const real_t density = (real_t)800.0;

        // A-frame tilt away from vertical (degrees)
        const real_t tiltDeg = (real_t)18.0;
        const real_t tiltRad = tiltDeg * (real_t)M_PI / (real_t)180.0;
        const real_t sinTilt = std::sin(tiltRad);
        const real_t cosTilt = std::cos(tiltRad);

        // Geometry: upright cards start flat, then rotate about Y by 90° ± tilt so their planes are parallel to yz
        // Feet are separated along x; all cards share the same y (yOffset).
        const real_t baseSpan = (real_t)4.0 * halfExtents.x() * std::tan(tiltRad) + (real_t)2.0 * halfExtents.z(); // small extra to avoid mid-height intersection
        const real_t centerHalfSpan = baseSpan * (real_t)0.5;

        // Layout parameters for a modest house (modifiable here)
        const int baseFrames = 10;    // A-frames on the bottom layer (placed along x)
        const int layers = baseFrames - 1;        // number of stacked layers of A-frames
        const real_t frameSpacing = baseSpan; //(real_t)0.085;  // spacing between adjacent frame centers along x
        const real_t roofGap = (real_t)0.00015;     // clearance above a ridge before placing a roof card
        const real_t layerGap = (real_t)0.00025;    // clearance between a roof top and the feet of the next layer
        const real_t layerSettleGap = (real_t)0.00005; // extra breathing room between layers so stacks can settle
        const real_t xOffset = (real_t)0.0;

        // Roof tilt: introduce a tiny pitch so the two supported ends differ by one card thickness.
        const real_t roofTilt = std::atan(halfExtents.z() / (baseSpan * 0.5));

        PhysicsSystem::CubeShape cardShape{halfExtents};

        auto placeLeaningCard = [&](real_t footX, real_t xCenter, real_t baseZ, real_t angleSign) {
            const real_t yaw = (real_t)M_PI_2 + angleSign * tiltRad; // rotate about Y: 90° ± tilt
            const real_t sinA = std::sin(yaw);
            const real_t cosA = std::cos(yaw);
            // vertical projection from long axis (x) and thickness (z)
            const real_t vertHalf = std::abs(cosA) * halfExtents.z() + std::abs(sinA) * halfExtents.x();
            Vector3r c(xOffset + xCenter + footX, (real_t)0.0, baseZ + vertHalf);
            Quaternion4r q = Quaternion4r(Eigen::AngleAxis<real_t>(yaw, Vector3r::UnitY()));
            PhysicsSystem::RigidState state; state.position = c; state.orientation = q; state.linearVelocity = Vector3r::Zero(); state.angularVelocity = Vector3r::Zero();
            PhysicsSystem::RigidProps props = PhysicsSystem::RigidProps::withDensity(density);
            sys.addRigidBody(cardShape, state, props);
            return baseZ + (real_t)2.0 * vertHalf;
        };

        auto placeAFrame = [&](real_t x, real_t baseZ) {
            const real_t leftTop = placeLeaningCard(-centerHalfSpan * 0.5, x, baseZ, (real_t)1.0);
            const real_t rightTop = placeLeaningCard(centerHalfSpan * 0.5, x, baseZ, (real_t)-1.0);
            return std::max(leftTop, rightTop);
        };

        auto placeRoof = [&](real_t x, real_t supportZ) {
            PhysicsSystem::RigidState state;
            state.position = Vector3r(xOffset + x, (real_t)0.0, supportZ + halfExtents.z() + roofGap);
            state.orientation = Quaternion4r(Eigen::AngleAxis<real_t>(roofTilt, Vector3r::UnitY()));
            state.linearVelocity = Vector3r::Zero();
            state.angularVelocity = Vector3r::Zero();
            PhysicsSystem::RigidProps props = PhysicsSystem::RigidProps::withDensity(density);
            sys.addRigidBody(cardShape, state, props);
        };

        real_t layerBaseZ = (real_t)0.0;
        int framesThisLayer = baseFrames;

        for (int layer = 0; layer < layers && framesThisLayer > 0; ++layer) {
            std::vector<real_t> xPositions(framesThisLayer);
            for (int i = 0; i < framesThisLayer; ++i) {
                xPositions[i] = ((real_t)i - ((real_t)framesThisLayer - (real_t)1.0) * (real_t)0.5) * frameSpacing;
            }

            real_t topZ = layerBaseZ;
            for (const real_t x : xPositions) {
                topZ = std::max(topZ, placeAFrame(x, layerBaseZ));
            }

            if (framesThisLayer > 1) {
                for (int j = 0; j + 1 < framesThisLayer; ++j) {
                    const real_t xMid = (xPositions[j] + xPositions[j + 1]) * (real_t)0.5;
                    placeRoof(xMid, topZ);
                }
                layerBaseZ = topZ + (real_t)2.0 * halfExtents.z() + roofGap + layerGap + layerSettleGap;
                --framesThisLayer;
            } else {
                break; // top ridge reached; no roof to support another layer
            }
        }

        // Ball rolling towards the house to knock it down
        const real_t ballRadius = (real_t)0.075;
        PhysicsSystem::SphereShape ballShape{ballRadius};
        PhysicsSystem::RigidState ballState;
        ballState.position = Vector3r(0.0, -1.0, ballRadius + 0.01);
        ballState.orientation = Quaternion4r::Identity();
        ballState.linearVelocity = Vector3r(0.0, 4.0, 0.0);
        ballState.angularVelocity = Vector3r::Zero();
        PhysicsSystem::RigidProps ballProps = PhysicsSystem::RigidProps::withDensity((real_t)1000.0);
        sys.addRigidBody(ballShape, ballState, ballProps);
    }
};
