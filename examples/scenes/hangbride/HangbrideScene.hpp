#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <vector>
#include "physics/constraints.hpp"

using namespace cardillo;

// HangbrideScene: Two opposing cliff blocks, two tripods per cliff (total 4),
// and two rope spans connecting tripod apex anchors across the gap.
class HangbrideScene : public SceneBase {
public:
    HangbrideScene() = default;
    ~HangbrideScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;
        using namespace cardillo::physics;
        const auto& cfg = sys.config();

        // Basic scales
        const real_t cliffHeight = (real_t)0.6;
        const real_t cliffThickness = (real_t)0.6;   // extent in x around its center
        const real_t cliffWidthY = (real_t)1.6;      // extent in y (half-width)
    const real_t gap = (real_t)7.2;              // gap between inner cliff faces (bridge length doubled again)
        const real_t cliffTopZ = (real_t)0.0;        // top surface at z=0
        const real_t cliffCenterZ = cliffTopZ - cliffHeight * (real_t)0.5;

        // Left and right cliff centers along X
        const real_t leftX  = -(gap * (real_t)0.5 + cliffThickness * (real_t)0.5);
        const real_t rightX =  (gap * (real_t)0.5 + cliffThickness * (real_t)0.5);

        // Create cliff blocks as static obstacle cubes
        auto makeCube = [&](const Vector3r& center, const Vector3r& halfExt, const Quaternion4r& q){
            PhysicsSystem::Cube c; c.center = center; c.halfExtents = halfExt; c.q = q; return c; };

        PhysicsSystem::Cube leftCliff = makeCube(Vector3r(leftX, 0.0, cliffCenterZ),
                                                 Vector3r(cliffThickness*(real_t)0.5, cliffWidthY*(real_t)0.5, cliffHeight*(real_t)0.5),
                                                 Quaternion4r::Identity());
        PhysicsSystem::Cube rightCliff = makeCube(Vector3r(rightX, 0.0, cliffCenterZ),
                                                  Vector3r(cliffThickness*(real_t)0.5, cliffWidthY*(real_t)0.5, cliffHeight*(real_t)0.5),
                                                  Quaternion4r::Identity());
        entt::entity eLeft  = entt::entity(static_cast<uint32_t>(sys.addObstacleBody(leftCliff)));
        entt::entity eRight = entt::entity(static_cast<uint32_t>(sys.addObstacleBody(rightCliff)));
        (void)eLeft; (void)eRight; // not needed further

        // Tripod parameters
        const real_t apexHeight = (real_t)0.6;  // height above top surface
        const real_t legFootRadius = (real_t)0.25; // spread around apex in x/y (unused)
        const real_t legSteeperInward = (real_t)0.22; // extra inward shift for cliff-facing leg to avoid overhang (unused)
        const real_t legK = (real_t)200000.0;    // stiff tripod springs
        const real_t legD = (real_t)2.0;        // tripod damping

        // Tripod anchor positions (apex) on left/right, two per side along Y
        const real_t apexInward = (real_t)0.05; // move slightly inward onto plateau from the inner edge
        std::vector<Vector3r> apexL = {
            Vector3r(leftX + cliffThickness*(real_t)0.5 - apexInward, -(real_t)0.5, cliffTopZ + apexHeight),
            Vector3r(leftX + cliffThickness*(real_t)0.5 - apexInward,  (real_t)0.5, cliffTopZ + apexHeight)
        };
        std::vector<Vector3r> apexR = {
            Vector3r(rightX - cliffThickness*(real_t)0.5 + apexInward, -(real_t)0.5, cliffTopZ + apexHeight),
            Vector3r(rightX - cliffThickness*(real_t)0.5 + apexInward,  (real_t)0.5, cliffTopZ + apexHeight)
        };

    // Build two tripods per side (left and right) as spring tripods
        // Returns apex point-mass entity for rope attachment
        auto buildTripod = [&](const Vector3r& apex, bool leftSide){
            // Create apex point mass (no rigid cube at top)
            auto createPointMass = [&](real_t mass, const Vector3r& p, real_t radius){
                auto& reg = sys.ecs();
                auto e = reg.create();
                reg.emplace<PhysicsSystem::C_PhysicsObject>(e);
                reg.emplace<PhysicsSystem::C_PointMassTag>(e);
                reg.emplace<PhysicsSystem::C_Collidable>(e);
                reg.emplace<PhysicsSystem::C_VisualObject>(e);
                reg.emplace<PhysicsSystem::C_PointVisualTag>(e);
                reg.emplace<PhysicsSystem::C_Mass>(e, PhysicsSystem::C_Mass{mass});
                reg.emplace<PhysicsSystem::C_Position3>(e, PhysicsSystem::C_Position3{p});
                reg.emplace<PhysicsSystem::C_LinearVelocity3>(e, PhysicsSystem::C_LinearVelocity3{Vector3r::Zero()});
                reg.emplace<PhysicsSystem::C_Radius>(e, PhysicsSystem::C_Radius{radius});
                // Friction component optional; omit to avoid registry version differences
                sys.markStructureDirty();
                return e;
            };

            entt::entity apexE = createPointMass((real_t)0.05, apex, (real_t)0.035);

            // Helper to add a 3D translational spring between apex point mass and obstacle cliff at local attachment point
            auto addLegSpring = [&](entt::entity apexEntity, entt::entity cliffEntity, const Vector3r& cliffLocal){
                sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), apexEntity, cliffEntity, Vector3r::Zero(), cliffLocal, legK, legD);
            };

            // Choose the relevant cliff and its local frame (identity rotation)
            const bool isLeft = leftSide;
            entt::entity cliffE = isLeft ? eLeft : eRight;

            // Compute local attachment points on cliff relative to the apex Y to avoid overlap between tripods:
            const real_t sx = (isLeft ? (real_t)+1.0 : (real_t)-1.0); // inner face is at local x = +halfX for left, -halfX for right
            const real_t halfX = cliffThickness * (real_t)0.5;
            const real_t halfY = cliffWidthY * (real_t)0.5;
            const real_t halfZ = cliffHeight * (real_t)0.5;

            // Center tripod under apex: two front legs at the inner top edge (±Y around apex.y), one back leg aligned in Y with apex and further inward in X.
            const real_t edgeMargin   = (real_t)0.005;  // attach on vertical inner face near top edge
            const real_t zDown        = (real_t)0.02;   // slightly below top edge on z
            const real_t frontSpreadY = (real_t)0.14;   // symmetric under apex (reduced to prevent overlap)
            const real_t backInward   = (real_t)0.26;   // stabilizing back leg further from edge

            // Apex Y in cliff-local frame (cliff centers are at y=0)
            const real_t yApexLocal = apex.y();
            auto clampY = [&](real_t y){ return std::min(std::max(y, -halfY + (real_t)0.02), halfY - (real_t)0.02); };

            // Local positions: front legs on inner vertical face near the top edge at y around apex.y; rear leg on top face inward at same y
            Vector3r pFrontPosLocal = Vector3r(sx * (halfX),  clampY(yApexLocal + frontSpreadY), halfZ - zDown);
            Vector3r pFrontNegLocal = Vector3r(sx * (halfX),  clampY(yApexLocal - frontSpreadY), halfZ - zDown);
            Vector3r pBackLocal     = Vector3r(sx * (halfX - backInward),  clampY(yApexLocal), halfZ);

            // Connect apex to these three points with stiff springs
            addLegSpring(apexE, cliffE, pFrontPosLocal);
            addLegSpring(apexE, cliffE, pFrontNegLocal);
            addLegSpring(apexE, cliffE, pBackLocal);

            return apexE; // use apex point mass as rope anchor
        };

        std::vector<entt::entity> anchorsL, anchorsR; anchorsL.reserve(2); anchorsR.reserve(2);
        anchorsL.push_back(buildTripod(apexL[0], true));
        anchorsL.push_back(buildTripod(apexL[1], true));
        anchorsR.push_back(buildTripod(apexR[0], false));
        anchorsR.push_back(buildTripod(apexR[1], false));

        // Helper: create a rope as a chain of point masses with springs; attach to two anchor entities
        auto addRope = [&](entt::entity eA, const Vector3r& pA, entt::entity eB, const Vector3r& pB,
                           int segments, real_t nodeMass, real_t k, real_t d){
            std::vector<entt::entity> nodes; nodes.reserve((size_t)segments);
            // Uniform initial positions along straight line
            for (int i = 0; i < segments; ++i) {
                real_t t = (segments == 1) ? (real_t)0.5 : (real_t)i / (real_t)(segments - 1);
                Vector3r p = (real_t)1.0 * ((real_t)1.0 - t) * pA + t * pB;
                (void)sys.addPointMass(nodeMass, p, Vector3r::Zero(), (real_t)0.03); // sphere radius
                // Retrieve the last created entity (point mass). addPointMass returns an index_t, not entity, so we can't directly get it here.
                // Instead, find the last physics object without C_BodyIndex yet would be brittle; better to create entities manually.
            }
        };

        // Since addPointMass returns index_t, create nodes explicitly to keep entity handles
        auto createPointMass = [&](real_t mass, const Vector3r& p, real_t radius){
            auto& reg = sys.ecs();
            auto e = reg.create();
            reg.emplace<PhysicsSystem::C_PhysicsObject>(e);
            reg.emplace<PhysicsSystem::C_PointMassTag>(e);
            reg.emplace<PhysicsSystem::C_Collidable>(e);
            reg.emplace<PhysicsSystem::C_VisualObject>(e);
            reg.emplace<PhysicsSystem::C_PointVisualTag>(e);
            reg.emplace<PhysicsSystem::C_Mass>(e, PhysicsSystem::C_Mass{mass});
            reg.emplace<PhysicsSystem::C_Position3>(e, PhysicsSystem::C_Position3{p});
            reg.emplace<PhysicsSystem::C_LinearVelocity3>(e, PhysicsSystem::C_LinearVelocity3{Vector3r::Zero()});
            reg.emplace<PhysicsSystem::C_Radius>(e, PhysicsSystem::C_Radius{radius});
            // Friction component optional; omit to avoid registry version differences
            sys.markStructureDirty();
            return e;
        };

        auto addRopeChain = [&](entt::entity eA, const Vector3r& pA, entt::entity eB, const Vector3r& pB,
                                int segments, real_t nodeMass, real_t k, real_t d) -> std::vector<entt::entity> {
            std::vector<entt::entity> nodes; nodes.reserve((size_t)segments);
            for (int i = 0; i < segments; ++i) {
                real_t t = (segments == 1) ? (real_t)0.5 : (real_t)i / (real_t)(segments - 1);
                Vector3r p = (real_t)1.0 * ((real_t)1.0 - t) * pA + t * pB;
                nodes.push_back(createPointMass(nodeMass, p, (real_t)0.03));
            }
            // Spring helper between two entities (point masses or anchor obstacles)
            auto addSpring = [&](entt::entity A, entt::entity B, const Vector3r& rA, const Vector3r& rB, real_t kmul){
                sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), A, B, rA, rB, k * kmul, d);
            };
            // Attach ends
            addSpring(nodes.front(), eA, Vector3r::Zero(), Vector3r::Zero(), (real_t)1.0);
            addSpring(nodes.back(),  eB, Vector3r::Zero(), Vector3r::Zero(), (real_t)1.0);
            // Chain segments
            for (int i = 0; i + 1 < segments; ++i) {
                addSpring(nodes[i], nodes[i+1], Vector3r::Zero(), Vector3r::Zero(), (real_t)1.0);
            }
            return nodes;
        };

    // Build two ropes across corresponding anchors
    const int segments = 120; // segment count for rope discretization
    const real_t ropeNodeMass = (real_t)0.02;
    const real_t ropeK = (real_t)3000000.0; // make top rope stiffer
    const real_t ropeD = (real_t)1.0;

        // Use the known apex positions for rope endpoints
        // Build ropes and keep their node entities for deck attachments
        std::vector<entt::entity> rope0 = addRopeChain(anchorsL[0], apexL[0], anchorsR[0], apexR[0], segments, ropeNodeMass, ropeK, ropeD);
        std::vector<entt::entity> rope1 = addRopeChain(anchorsL[1], apexL[1], anchorsR[1], apexR[1], segments, ropeNodeMass, ropeK, ropeD);

        // Helper to add a 3D translational spring
        auto add3DSpring = [&](entt::entity A, entt::entity B, const Vector3r& rA, const Vector3r& rB, real_t k, real_t d){
            sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), A, B, rA, rB, k, d);
        };

    // Add floor boards (planks) hanging from the two ropes
    const int marginSegments = 2;                  // symmetric margin of rope nodes at both ends
    const int iStart = marginSegments;
    const int iEnd   = (segments - 1) - marginSegments;
    const int stride = 4;                          // increase spacing between planks to create a gap
        const real_t plankMass = (real_t)0.35;
        // Long along span (x), wide across ropes (y), thinner boards (z)
        Vector3r plankHalf = Vector3r((real_t)0.10, (real_t)0.46, (real_t)0.02);
        const real_t plankHangK = ropeK, plankHangD = ropeD;
        const real_t plankLinkK = ropeK, plankLinkD = ropeD;

        // Precompute straight-line param for initial placement
        auto posAlong = [&](const Vector3r& A, const Vector3r& B, real_t t){ return ((real_t)1 - t) * A + t * B; };

        std::vector<entt::entity> planks; planks.reserve((size_t)std::max(0, iEnd - iStart + 2));
        // Build the list of rope indices for planks with equalized end margins within one stride
        std::vector<int> plankIndices;
        plankIndices.reserve((size_t)std::max(1, (iEnd - iStart) / stride + 2));
        {
            const int remainder = (iEnd - iStart) % stride;
            const int leftPad  = remainder / 2;
            const int rightPad = remainder - leftPad;
            const int iStartAdj = iStart + leftPad;
            const int iEndAdj   = iEnd   - rightPad;
            for (int i = iStartAdj; i <= iEndAdj; i += stride) plankIndices.push_back(i);
            if (plankIndices.empty()) plankIndices.push_back((iStart + iEnd) / 2); // fallback: center one
        }

        // Precompute inner face x-positions and a small clearance/lift to avoid initial intersections
        const real_t leftInnerX  = leftX  + cliffThickness * (real_t)0.5;
        const real_t rightInnerX = rightX - cliffThickness * (real_t)0.5;
        const real_t edgeClearX  = (real_t)0.01; // 1 cm clearance beyond plank half-length
        const real_t liftZ       = (real_t)0.004; // 4 mm above top surface to avoid initial clipping

        for (int i : plankIndices) {
            real_t t = (real_t)i / (real_t)(segments - 1);
            // Initial centers along the straight lines between apexes, but positioned near the height of the cliff floor (bottom)
            Vector3r p0 = posAlong(apexL[0], apexR[0], t);
            Vector3r p1 = posAlong(apexL[1], apexR[1], t);
            Vector3r center = (p0 + p1) * (real_t)0.5;
            // Clamp along X so the board does not intersect the cliff inner faces at either end
            const real_t minX = leftInnerX  + plankHalf.x() + edgeClearX;
            const real_t maxX = rightInnerX - plankHalf.x() - edgeClearX;
            if (minX <= maxX) {
                center.x() = std::min(std::max(center.x(), minX), maxX);
            }
            // Place so the top surface sits slightly above the cliffs to avoid initial penetration
            center.z() = (cliffTopZ - plankHalf.z()) + liftZ;

            // Create plank as a rigid body cube (identity orientation)
            PhysicsSystem::Cube shape; shape.center = center; shape.halfExtents = plankHalf; shape.q = Quaternion4r::Identity();
            entt::entity plank = sys.addRigidBody(plankMass, center, Quaternion4r::Identity(), Vector3r::Zero(), Vector3r::Zero(), shape);
            planks.push_back(plank);

            // Attach plank to rope nodes at its top-left and top-right surfaces
            add3DSpring(plank, rope0[i], Vector3r(0, -plankHalf.y(),  plankHalf.z()), Vector3r::Zero(), plankHangK, plankHangD);
            add3DSpring(plank, rope1[i], Vector3r(0,  plankHalf.y(),  plankHalf.z()), Vector3r::Zero(), plankHangK, plankHangD);
        }

        // Connect adjacent planks at their left and right surfaces to form a continuous deck
        // Replace single stiff springs between neighboring planks by a short rope (series of springs)
        // so the connection behaves more like a rope (limits extension more than contraction).
        auto addShortRopeBetween = [&](entt::entity A, entt::entity B,
                                       const Vector3r& rA, const Vector3r& rB,
                                       int segmentsRope, real_t kRope, real_t dRope) {
            // Compute approximate world positions for A and B attachment points.
            // Use position component if available; if not, fall back to origin.
            auto& reg = sys.ecs();
            Vector3r posA = Vector3r::Zero();
            Vector3r posB = Vector3r::Zero();
            if (auto pa = reg.try_get<PhysicsSystem::C_Position3>(A)) posA = pa->value + rA;
            else posA = rA;
            if (auto pb = reg.try_get<PhysicsSystem::C_Position3>(B)) posB = pb->value + rB;
            else posB = rB;

            // Create rope nodes (point masses) between posA and posB
            std::vector<entt::entity> nodes; nodes.reserve((size_t)segmentsRope);
            for (int s = 0; s < segmentsRope; ++s) {
                real_t t = (segmentsRope == 1) ? (real_t)0.5 : (real_t)s / (real_t)(segmentsRope - 1);
                Vector3r p = ((real_t)1 - t) * posA + t * posB;
                nodes.push_back(createPointMass((real_t)0.02, p, (real_t)0.02));
            }

            // helper to add translational springs between two entities
            auto addLink = [&](entt::entity X, entt::entity Y, const Vector3r& rx, const Vector3r& ry){
                sys.addConstraint<cardillo::physics::LinearDistanceConstraint>(sys.ecs(), X, Y, rx, ry, kRope, dRope);
            };

            // Attach ends
            addLink(A, nodes.front(), rA, Vector3r::Zero());
            addLink(nodes.back(), B, Vector3r::Zero(), rB);

            // Chain internal nodes
            for (int s = 0; s + 1 < segmentsRope; ++s) {
                addLink(nodes[s], nodes[s+1], Vector3r::Zero(), Vector3r::Zero());
            }
        };

        for (size_t j = 0; j + 1 < planks.size(); ++j) {
            entt::entity A = planks[j];
            entt::entity B = planks[j+1];
            // Left edges (toward rope0): connect end faces (x surfaces) at top-left corners using a 3-segment rope
            addShortRopeBetween(A, B,
                               Vector3r( plankHalf.x(), -plankHalf.y(),  plankHalf.z()),
                               Vector3r(-plankHalf.x(), -plankHalf.y(),  plankHalf.z()),
                               3, plankLinkK, plankLinkD);
            // Right edges (toward rope1): connect end faces (x surfaces) at top-right corners using a 3-segment rope
            addShortRopeBetween(A, B,
                               Vector3r( plankHalf.x(),  plankHalf.y(),  plankHalf.z()),
                               Vector3r(-plankHalf.x(),  plankHalf.y(),  plankHalf.z()),
                               3, plankLinkK, plankLinkD);
        }

        // Connect the deck's outer ends to the cliff inner edges for added stability
        if (!planks.empty()) {
            const real_t halfX = cliffThickness * (real_t)0.5;
            const real_t halfY = cliffWidthY * (real_t)0.5; (void)halfY; // not strictly needed but informative
            const real_t halfZ = cliffHeight * (real_t)0.5;
            const real_t edgeZDown = (real_t)0.02; // attach slightly below the top edge on the vertical face
            const real_t kEdge = (real_t)4000.0;
            const real_t dEdge = (real_t)0.35;

            entt::entity firstPlank = planks.front();
            entt::entity lastPlank  = planks.back();

            // Left cliff: inner face is at local x = +halfX
            add3DSpring(firstPlank, eLeft,
                        Vector3r(-plankHalf.x(), -plankHalf.y(),  plankHalf.z()),
                        Vector3r(+halfX,         -plankHalf.y(),  halfZ - edgeZDown),
                        kEdge, dEdge);
            add3DSpring(firstPlank, eLeft,
                        Vector3r(-plankHalf.x(),  plankHalf.y(),  plankHalf.z()),
                        Vector3r(+halfX,          plankHalf.y(),  halfZ - edgeZDown),
                        kEdge, dEdge);

            // Right cliff: inner face is at local x = -halfX
            add3DSpring(lastPlank, eRight,
                        Vector3r(+plankHalf.x(), -plankHalf.y(),  plankHalf.z()),
                        Vector3r(-halfX,         -plankHalf.y(),  halfZ - edgeZDown),
                        kEdge, dEdge);
            add3DSpring(lastPlank, eRight,
                        Vector3r(+plankHalf.x(),  plankHalf.y(),  plankHalf.z()),
                        Vector3r(-halfX,          plankHalf.y(),  halfZ - edgeZDown),
                        kEdge, dEdge);
        }
    }
};
