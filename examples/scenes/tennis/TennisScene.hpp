#pragma once

#include "../SceneBase.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "misc/spline.hpp"

using namespace cardillo;

class TennisScene : public SceneBase {
   public:
    const char* sceneName() const override { return "tennis"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        engine.setGravity(Vector3r(-9.81, 0, 0));

        // ── Frame geometry ────────────────────────────────────────────────────
        const real_t frameHalfX = (real_t)0.010;
        const real_t innerHalfY = (real_t)0.130;
        const real_t innerHalfZ = (real_t)0.130;
        const real_t frameBar = (real_t)0.008;

        const Vector3r topBottomHalf(frameHalfX, innerHalfY + frameBar, frameBar);
        const Vector3r sideHalf(frameHalfX, frameBar, innerHalfZ);

        const Vector3r topPos((real_t)0.0, (real_t)0.0, innerHalfZ + frameBar);
        const Vector3r bottomPos((real_t)0.0, (real_t)0.0, -(innerHalfZ + frameBar));
        const Vector3r leftPos((real_t)0.0, -(innerHalfY + frameBar), (real_t)0.0);
        const Vector3r rightPos((real_t)0.0, (innerHalfY + frameBar), (real_t)0.0);

        const physics::RigidProps frameProps = physics::RigidProps::withDensity((real_t)650.0);
        const entt::entity frameTop = engine.addRigidBody(physics::CubeShape(topBottomHalf), physics::RigidState(topPos), frameProps);
        const entt::entity frameBottom = engine.addRigidBody(physics::CubeShape(topBottomHalf), physics::RigidState(bottomPos), frameProps);
        const entt::entity frameLeft = engine.addRigidBody(physics::CubeShape(sideHalf), physics::RigidState(leftPos), frameProps);
        const entt::entity frameRight = engine.addRigidBody(physics::CubeShape(sideHalf), physics::RigidState(rightPos), frameProps);

        const Vector3r handleHalf((real_t)0.012, (real_t)0.018, (real_t)0.140);
        const Vector3r handlePos((real_t)0.0, (real_t)0.0, -(innerHalfZ + (real_t)2.0 * frameBar + handleHalf.z()));
        const entt::entity handle = engine.addRigidBody(physics::CubeShape(handleHalf), physics::RigidState(handlePos), physics::RigidProps::withDensity((real_t)520.0));

        engine.addRigidConstraint(frameTop, handle);
        engine.addRigidConstraint(frameBottom, handle);
        engine.addRigidConstraint(frameLeft, handle);
        engine.addRigidConstraint(frameRight, handle);

        engine.disableCollisionBetween(frameTop, frameBottom);
        engine.disableCollisionBetween(frameTop, frameLeft);
        engine.disableCollisionBetween(frameTop, frameRight);
        engine.disableCollisionBetween(frameBottom, frameLeft);
        engine.disableCollisionBetween(frameBottom, frameRight);
        engine.disableCollisionBetween(frameLeft, frameRight);
        engine.disableCollisionBetween(frameTop, handle);
        engine.disableCollisionBetween(frameBottom, handle);
        engine.disableCollisionBetween(frameLeft, handle);
        engine.disableCollisionBetween(frameRight, handle);

        const Vector3r Ktrans((real_t)250000.0, (real_t)250000.0, (real_t)250000.0);
        const Vector3r Dtrans((real_t)0.0, (real_t)0.0, (real_t)0.0);
        const Vector3r Krot((real_t)45000.0, (real_t)45000.0, (real_t)45000.0);
        const Vector3r Drot((real_t)0.0, (real_t)0.0, (real_t)0.0);
        engine.addTranslationRotationConstraint(handle, entt::null, physics::JointFrame(handle), Ktrans, Dtrans, Krot, Drot);

        // ── String grid ───────────────────────────────────────────────────────
        const int totalMains = 18;
        const int totalCrosses = 21;
        const int trimEdge = 1;

        const real_t stringRadius = (real_t)0.000625;  // 1.25 mm gauge

        // zigAmp: each layer offsets by exactly 1× radius so capsule surfaces
        // just touch (total separation = 2r = 1 diameter). No extra factor needed.
        const real_t zigAmp = stringRadius * 0.5;

        const real_t yMin = -innerHalfY;
        const real_t yMax = innerHalfY;
        const real_t zMin = -innerHalfZ;
        const real_t zMax = innerHalfZ;

        const auto lerp = [](real_t a, real_t b, real_t t) { return a + (b - a) * t; };

        // Crossing positions — shared between both layers.
        std::vector<real_t> mainY(totalMains);
        for (int col = 0; col < totalMains; ++col) mainY[col] = lerp(yMin, yMax, (real_t)col / (real_t)(totalMains - 1));

        std::vector<real_t> crossZ(totalCrosses);
        for (int row = 0; row < totalCrosses; ++row) crossZ[row] = lerp(zMin, zMax, (real_t)row / (real_t)(totalCrosses - 1));

        const physics::BeamCrossSection stringSection(stringRadius, stringRadius, physics::BeamBodyType::Capsule);
        physics::BeamSpringParams stringSprings =
            physics::BeamSpringParams::fromMaterial((real_t)5.0e10, (real_t)0.35, (real_t)1.0, (real_t)1.0, (real_t)1.0, (real_t)1.0, (real_t)1.0, (real_t)0.0005);
        stringSprings.gamma0 = Vector3r::Zero();
        stringSprings.gamma0->x() = (real_t)-1.0e-5;
        const physics::RigidProps stringProps = physics::RigidProps::withDensity((real_t)1200.0);
        const physics::RigidState stringStateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero());

        auto buildPolylineBeam = [&](const std::vector<Vector3r>& pts, entt::entity attachA, entt::entity attachB, const std::string& tag) {
            std::vector<std::unique_ptr<LinearSpline>> ownedParts;
            std::vector<const SplinePattern*> parts;
            ownedParts.reserve(pts.size());
            parts.reserve(pts.size());
            for (size_t i = 1; i < pts.size(); ++i) {
                if ((pts[i] - pts[i - 1]).norm() <= (real_t)1e-10) continue;
                ownedParts.push_back(std::make_unique<LinearSpline>(pts[i - 1], pts[i]));
                parts.push_back(ownedParts.back().get());
            }
            if (parts.empty()) return;
            auto ends = engine.createBeams(parts, stringSection, stringSprings, stringStateDefaults, stringProps, (size_t)64);
            if (ends.first != entt::null) {
                engine.addRigidConstraint(ends.first, attachA);
                engine.disableCollisionBetween(ends.first, attachA);
            }
            if (ends.second != entt::null) {
                engine.addRigidConstraint(ends.second, attachB);
                engine.disableCollisionBetween(ends.second, attachB);
            }
        };

        // ── Main strings: run along Z, fixed Y ───────────────────────────────

        const real_t halfGapZ = 0;
        const real_t halfGapY = 0;

        for (int col = 0; col < totalMains; ++col) {
            if (col < trimEdge || col >= totalMains - trimEdge) continue;

            const real_t y = mainY[col];

            std::vector<Vector3r> pts;
            pts.reserve(totalCrosses + 4);

            // Anchor + approach: both at x=0, horizontally level with the frame face
            pts.push_back(Vector3r((real_t)0.0, y, crossZ[0] - halfGapZ));

            // One node per crossing — direct alternating ±zigAmp, no midpoints
            for (int row = 0; row < totalCrosses; ++row) {
                const int sign = ((col + row) % 2 == 0) ? 1 : -1;
                pts.push_back(Vector3r(zigAmp * (real_t)sign, y, crossZ[row]));
            }

            // Runout + anchor: mirror of approach
            pts.push_back(Vector3r((real_t)0.0, y, crossZ[totalCrosses - 1] + halfGapZ));

            buildPolylineBeam(pts, frameBottom, frameTop, "string_v_" + std::to_string(col));
        }

        // ── Cross strings: run along Y, fixed Z ──────────────────────────────
        // Mirror of above. Parity +1 so crosses are always on the opposite X
        // side from mains at each crossing → genuine over/under weave.

        for (int row = 0; row < totalCrosses; ++row) {
            if (row < trimEdge || row >= totalCrosses - trimEdge) continue;

            const real_t z = crossZ[row];

            std::vector<Vector3r> pts;
            pts.reserve(totalMains + 4);

            pts.push_back(Vector3r((real_t)0.0, mainY[0] - halfGapY, z));

            for (int col = 0; col < totalMains; ++col) {
                const int sign = ((col + row + 1) % 2 == 0) ? 1 : -1;
                pts.push_back(Vector3r(zigAmp * (real_t)sign, mainY[col], z));
            }

            pts.push_back(Vector3r((real_t)0.0, mainY[totalMains - 1] + halfGapY, z));

            buildPolylineBeam(pts, frameLeft, frameRight, "string_h_" + std::to_string(row));
        }

        // ── Tennis ball ───────────────────────────────────────────────────────
        physics::RigidState ballState;
        ballState.position = Vector3r((real_t)-0.05, (real_t)0.0, (real_t)0.0);
        ballState.linearVelocity = Vector3r((real_t)7.0, (real_t)0.0, (real_t)0.0);
        // engine.addSoftBody("res/meshes/sphere.obj", (real_t)1.0e5, (real_t)0.05, ballState.position, ballState.orientation, ballState.linearVelocity, ballState.angularVelocity, (real_t)0.057,
        //                    (real_t)0.005);
        engine.addRigidBody(physics::SphereShape((real_t)0.033), ballState, physics::RigidProps::withDensity((real_t)800.0));
    }
};