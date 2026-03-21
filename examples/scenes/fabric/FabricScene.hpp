#pragma once

#include "../SceneBase.hpp"
#include "misc/spline.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <filesystem>
#include <limits>
#include <unordered_set>
#include <vector>
#include <unistd.h>

using namespace cardillo;

class FabricScene : public SceneBase {
public:
    const char* sceneName() const override { return "fabric"; }
    FabricScene() = default;
    ~FabricScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        auto& sys = engine.world();
        using namespace cardillo;
        using namespace cardillo::misc;

        sys.setGravity(Vector3r(0, 0, -9.81));

        // Plane
        // engine.addStaticBody(physics::PlaneShape(Vector3r(0, 0, 1), Vector3r(0, 1, 0), (real_t)8, (real_t)8),
        //                   physics::RigidState(Vector3r(0, 0, 0)));

        // Cone
        // engine.addStaticBody(physics::ConeShape((real_t)0.01, (real_t)0.10),
        //             physics::RigidState(Vector3r(0, 0, 0)));

        // Sphere
        // engine.addStaticBody(physics::SphereShape((real_t)0.05),
        //             physics::RigidState(Vector3r(0, 0, 0)));

        // Load BCC splines
        const std::filesystem::path bccRel = "res/bcc/openwork_trellis_pattern.bcc";  //"res/bcc/openwork_trellis_pattern.bcc" "res/bcc/flame_ribbing_pattern.bcc"
        std::vector<std::filesystem::path> candidates = {
            bccRel,
            std::filesystem::path("../") / bccRel,
            std::filesystem::path("../../") / bccRel,
            std::filesystem::path("../../../") / bccRel,
        };
        const std::filesystem::path here = std::filesystem::path(__FILE__).parent_path();
        const std::filesystem::path root = here.parent_path().parent_path().parent_path().parent_path();
        candidates.push_back(root / bccRel);
        {
            std::array<char, PATH_MAX> buf{};
            const ssize_t len = ::readlink("/proc/self/exe", buf.data(), buf.size() - 1);
            if (len > 0) {
                buf[(size_t)len] = '\0';
                std::filesystem::path exeDir = std::filesystem::path(buf.data()).parent_path();
                std::filesystem::path exeRoot = exeDir.parent_path().parent_path();
                candidates.push_back(exeRoot / bccRel);
            }
        }

        std::filesystem::path bccPath;
        for (const auto& c : candidates) {
            if (std::filesystem::exists(c)) { bccPath = c; break; }
        }

        const real_t bccScale = (real_t)0.01;
        auto splines = loadSplinesFromBCC(bccPath.empty() ? bccRel.string() : bccPath.string(), bccScale);
        if (splines.empty()) {
            const std::string pathStr = bccPath.empty() ? bccRel.string() : bccPath.string();
            std::cerr << "[FabricScene] Failed to load BCC: " << pathStr << std::endl;
            return;
        }

        // Compute bounding box by sampling each spline
        Vector3r bbMin(std::numeric_limits<real_t>::infinity(),
                       std::numeric_limits<real_t>::infinity(),
                       std::numeric_limits<real_t>::infinity());
        Vector3r bbMax(-std::numeric_limits<real_t>::infinity(),
                       -std::numeric_limits<real_t>::infinity(),
                       -std::numeric_limits<real_t>::infinity());
        const int samplesPerSpline = 64;
        for (const auto& sp : splines) {
            if (!sp) continue;
            for (int i = 0; i <= samplesPerSpline; ++i) {
                real_t a = (real_t)i / (real_t)samplesPerSpline;
                Vector3r p = sp->sample(a).position;
                bbMin = bbMin.cwiseMin(p);
                bbMax = bbMax.cwiseMax(p);
            }
        }

        // Print bounding box size
        const Vector3r bbSize = bbMax - bbMin;
        std::cout << "[FabricScene] BCC bounding box size: " << bbSize.transpose() << std::endl;

        const Vector3r centerXY((bbMin.x() + bbMax.x()) * (real_t)0.5,
                    (bbMin.y() + bbMax.y()) * (real_t)0.5,
                    (real_t)0.0);
        const real_t lift = (real_t)0.07;
        const Vector3r offset = Vector3r(-centerXY.x(), -centerXY.y(), -bbMin.z() + lift);
        std::cout << "[FabricScene] Applying offset to splines: " << offset.transpose() << std::endl;

        // String material properties
        const real_t diameter = (real_t)0.00175; // m
        const real_t density  = (real_t)800.0; // kg/m^3 
        const real_t E        = (real_t)5e6;    // Pa
        const real_t nu       = (real_t)0.40;

        physics::BeamCrossSection section(diameter, diameter, physics::BeamBodyType::Capsule);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu, (real_t)1, (real_t)1, (real_t)0.05, (real_t)0.05, (real_t)0.05);
        // springs.setDampingFromFactor((real_t)0.02);

        physics::RigidProps props = physics::RigidProps::withDensity(density);
        physics::RigidState stateDefaults(offset, Vector3r(0, 0, 0));

        const float totalLength = std::accumulate(splines.begin(), splines.end(), 0.0f,
                                              [](float sum, const std::shared_ptr<SplinePattern>& sp) {
                                                  return sum + (sp ? sp->totalLength() : 0.0f);
                                              });

        const int resampleSegments = totalLength / (diameter * 1.05);
       

        std::cout << "[FabricScene] Total spline length: " << totalLength << " m, average segment length: "
                  << (totalLength / (real_t)resampleSegments) << " m" << std::endl;
        for (const auto& sp : splines) {
            if (!sp) continue;
            int segments = sp.get()->totalLength() / (totalLength / (real_t)resampleSegments);
            engine.createBeam(*sp.get(), section, springs, stateDefaults, props, segments);
            std::cout << "[FabricScene] Created beam from spline with " << segments << " segments." << std::endl;
        }

        // Close open loops manually by rigidly constraining ends to elements 6 steps inward for "res/bcc/openwork_trellis_pattern.bcc"
        if(bccRel == "res/bcc/openwork_trellis_pattern.bcc") {
            auto& reg = sys.ecs();
            auto view = reg.view<World::C_BeamElement>();
            std::unordered_set<entt::entity> visited;
            const size_t offsetSteps = 6;

            for (auto e : view) {
                if (visited.count(e)) continue;
                // Find chain start by walking prev links
                entt::entity start = e;
                bool isLoop = false;
                std::unordered_set<entt::entity> backSeen;
                while (reg.valid(start) && reg.any_of<World::C_BeamElement>(start)) {
                    const auto& be = reg.get<World::C_BeamElement>(start);
                    if (!be.prev.has_value()) break;
                    if (be.prev.value() == start || backSeen.count(be.prev.value())) { isLoop = true; break; }
                    backSeen.insert(start);
                    start = be.prev.value();
                }

                // If this is a loop (no end), skip
                if (isLoop) {
                    continue;
                }

                // Collect the chain
                std::vector<entt::entity> chain;
                entt::entity cur = start;
                while (reg.valid(cur) && reg.any_of<World::C_BeamElement>(cur)) {
                    if (visited.count(cur)) break;
                    visited.insert(cur);
                    chain.push_back(cur);
                    const auto& be = reg.get<World::C_BeamElement>(cur);
                    if (!be.next.has_value()) break;
                    if (be.next.value() == cur) break;
                    cur = be.next.value();
                }

                if (chain.size() <= offsetSteps) continue;

                const entt::entity end = chain.back();
                if (!reg.valid(end)) continue;

                // Constrain start to element 6 prior to end, and end to element 6 next from start
                const entt::entity startNext = chain[offsetSteps];
                const entt::entity endPrev = chain[chain.size() - 1 - offsetSteps];

                if (reg.valid(start) && reg.valid(endPrev)) {
                    sys.addRigidConstraint(start, endPrev);
                    sys.disableCollisionBetween(start, endPrev);
                }
                if (reg.valid(end) && reg.valid(startNext)) {
                    sys.addRigidConstraint(end, startNext);
                    sys.disableCollisionBetween(end, startNext);
                }
            }
        }

        // Compute world-space AABB after offset
        const Vector3r bbMinW = bbMin + offset;
        const Vector3r bbMaxW = bbMax + offset;

        // Create two super-heavy cubes at +/- X edges of the fabric AABB
        const real_t cubeThickness = std::max((real_t)0.005, bbSize.x() * (real_t)0.05);
        const Vector3r cubeHalfExtents(cubeThickness,
                                       bbSize.y() * (real_t)0.5,
                                       bbSize.z() * (real_t)0.5);
        const real_t centerY = (bbMinW.y() + bbMaxW.y()) * (real_t)0.5;
        const real_t centerZ = (bbMinW.z() + bbMaxW.z()) * (real_t)0.5;

        const Vector3r leftCubePos(bbMinW.x() - cubeHalfExtents.x() - 0.0025, centerY, centerZ);
        const Vector3r rightCubePos(bbMaxW.x() + cubeHalfExtents.x() + 0.0075, centerY, centerZ);

        m_leftCube = engine.addRigidBody(physics::CubeShape(cubeHalfExtents),
                                      physics::RigidState(leftCubePos),
                                      physics::RigidProps((real_t)1e10));
        m_rightCube = engine.addRigidBody(physics::CubeShape(cubeHalfExtents),
                                       physics::RigidState(rightCubePos),
                                       physics::RigidProps((real_t)1e10));

        // Select 10 beam elements on each +/-X edge, closest to targets on the inner face
        struct Candidate { entt::entity e; Vector3r p; };
        std::vector<Candidate> allCandidates;

        auto view = sys.ecs().view<World::C_BeamElement, World::C_Position3>();
        real_t minX = std::numeric_limits<real_t>::infinity();
        real_t maxX = -std::numeric_limits<real_t>::infinity();
        for (auto [e, be, pos] : view.each()) {
            (void)e; (void)be;
            const real_t x = pos.value.x();
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
        }
        for (auto [e, be, pos] : view.each()) {
            (void)be;
            allCandidates.push_back({e, pos.value});
        }

        auto pickClosestToTargets = [&](real_t xTarget, size_t count) {
            std::vector<entt::entity> selected;
            if (allCandidates.empty()) return selected;

            std::vector<bool> used(allCandidates.size(), false);
            const real_t yMin = bbMinW.y();
            const real_t yMax = bbMaxW.y();
            const real_t span = std::max((real_t)1e-6, yMax - yMin);
            const real_t zTarget = (bbMinW.z() + bbMaxW.z()) * (real_t)0.5;

            for (size_t i = 0; i < count; ++i) {
                const real_t yTarget = yMin + (span * ((real_t)i + (real_t)0.5) / (real_t)count);
                const Vector3r target(xTarget, yTarget, zTarget);
                real_t bestDist = std::numeric_limits<real_t>::infinity();
                size_t bestIdx = allCandidates.size();
                for (size_t j = 0; j < allCandidates.size(); ++j) {
                    if (used[j]) continue;
                    const real_t d = (allCandidates[j].p - target).squaredNorm();
                    if (d < bestDist) { bestDist = d; bestIdx = j; }
                }
                if (bestIdx < allCandidates.size()) {
                    used[bestIdx] = true;
                    selected.push_back(allCandidates[bestIdx].e);
                }
            }
            return selected;
        };

        const size_t attachCount = 10;
        auto leftAttach = pickClosestToTargets(minX, attachCount);
        auto rightAttach = pickClosestToTargets(maxX, attachCount);

        for (const auto& e : leftAttach) {
            sys.addRigidConstraint(m_leftCube, e);
            sys.disableCollisionBetween(m_leftCube, e);
        }
        for (const auto& e : rightAttach) {
            sys.addRigidConstraint(m_rightCube, e);
            sys.disableCollisionBetween(m_rightCube, e);
        }

        // Motion parameters based on AABB size
        m_moveDuration = (real_t)0.5;
        m_twistRamp = (real_t)0.5;
        m_moveSpeed = bbSize.x() * (real_t)0.15 / m_moveDuration; // toward each other
        m_twistOmega = (real_t)2.0; // rad/s around Z

        // Print total number of beam elements created
        size_t beamCount = 0;
        sys.ecs().view<World::C_BeamElement>().each([&](auto&) { ++beamCount; });
        // std::cout << "[FabricScene] Created " << beamCount << " beam elements." << std::endl;
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        auto& sys = engine.world();
        using namespace cardillo;

        auto smoothstep = [](real_t x) {
            x = std::clamp(x, (real_t)0, (real_t)1);
            return x * x * ((real_t)3 - (real_t)2 * x);
        };

        const bool haveCubes = (m_leftCube != entt::null) && (m_rightCube != entt::null);
        if (!haveCubes) return;

        // Cancel gravity
        sys.applyForce(m_leftCube,  -sys.gravity() * sys.getMass(m_leftCube).diagonal(), Vector3r::Zero());
        sys.applyForce(m_rightCube, -sys.gravity() * sys.getMass(m_rightCube).diagonal(), Vector3r::Zero());

        // First phase: move toward each other for 0.5 s with smooth accel/decel
        if (t <= m_moveDuration) {
            const real_t s = t / m_moveDuration;
            const real_t v = m_moveSpeed * std::sin((real_t)M_PI * s);
            sys.setLinearVelocity(m_leftCube,  Vector3r(+v, 0, 0));
            sys.setLinearVelocity(m_rightCube, Vector3r(-v, 0, 0));
            sys.setAngularVelocity(m_leftCube,  Vector3r::Zero());
            sys.setAngularVelocity(m_rightCube, Vector3r::Zero());
            return;
        }

        // Second phase: twist around Z in opposite directions, ramped up smoothly
        const real_t sTwist = (t - m_moveDuration) / std::max((real_t)1e-6, m_twistRamp);
        const real_t omega = m_twistOmega * smoothstep(sTwist);
        sys.setLinearVelocity(m_leftCube,  Vector3r::Zero());
        sys.setLinearVelocity(m_rightCube, Vector3r::Zero());
        sys.setAngularVelocity(m_leftCube,  Vector3r(0, +omega, 0));
        sys.setAngularVelocity(m_rightCube, Vector3r(0, -omega, 0));
    }

private:
    entt::entity m_leftCube{entt::null};
    entt::entity m_rightCube{entt::null};
    real_t m_moveDuration{(real_t)0.5};
    real_t m_twistRamp{(real_t)0.5};
    real_t m_moveSpeed{(real_t)0.1};
    real_t m_twistOmega{(real_t)2.0};
};
