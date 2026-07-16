#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <optional>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

using namespace cardillo;

// Strandbeest single-leg scene: builds a static bar with anchors p1, p2 and a single
// Theo-Jansen inspired leg constructed from planar link lengths.
class StrandbeestScene : public SceneBase {
public:
    const char* sceneName() const override { return "strandbeest"; }
    StrandbeestScene() = default;
    ~StrandbeestScene() override = default;
    
    inline static const real_t stiffness = 1e6;
    inline static const real_t damping   = 1e6;
    inline static const Vector3r K_TRANS = Vector3r::Constant(stiffness);
    inline static const Vector3r D_TRANS = Vector3r::Constant(damping);
    inline static const Vector3r R_ROT = Vector3r((real_t)0.0, damping, damping);
    inline static const Vector3r D_ROT = Vector3r((real_t)1e-3, damping, damping);

    struct LegParams {
        real_t dx{0.380};
        real_t dz{0.078};
    };

    // The holy numbers from Theo-Jansen himself
    struct Lengths {
        real_t b{0.415};
        real_t c{0.393};
        real_t d{0.401};
        real_t e{0.558};
        real_t f{0.394};
        real_t g{0.367};
        real_t h{0.657};
        real_t i{0.490};
        real_t j{0.500};
        real_t k{0.619};
        real_t m{0.150};
    };

    using SegDef = std::pair<std::string, std::string>;
    const std::unordered_map<std::string, SegDef> kSegments {
        {"frame", {"p1",    "p2"   }},
        {"m",     {"p2",    "mjk"  }},
        {"k",     {"gcki",   "mjk"  }},
        {"j",     {"ebj",    "mjk"  }},
        {"b",     {"p1",    "ebj"  }},
        {"e",     {"edf",    "ebj"  }},
        {"f",     {"edf",   "fgh"  }},
        {"d",     {"edf",    "p1"  }},
        {"c",     {"p1",   "gcki" }},
        {"g",     {"gcki",  "fgh" }},
        {"h",     {"fgh",    "hi"  }},
        {"i",     {"gcki",   "hi"   }},
    };

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        const real_t slopeDeg = (real_t)3.0;
        const real_t slopeRad = slopeDeg * (real_t)M_PI / (real_t)180.0;
        // Tilt gravity instead of the ground plane: downhill +x
        engine.setGravity(Vector3r(-9.81 * std::sin(slopeRad), 0.0, -9.81 * std::cos(slopeRad)));

        // Configurable params
        LegParams params;
        Lengths L;
        const real_t baseAlphaDeg = 20.0; // crank angle for bar m (degrees)
        const real_t walkerLift = (real_t)0.85 + 0.01; // raise walker to clear ground
        const Vector3r p2_world(params.dx, 0.0, params.dz + walkerLift);
        const Vector3r planeNormal = Vector3r::UnitY(); // keep walker upright; leg plane normal
        const Vector3r xAxis = Vector3r::UnitX();
        const int nLegPairs = 8;
        // Ground cube with top surface at z=0 (replaces plane)
        physics::CubeShape groundShape{Vector3r((real_t)20.0, (real_t)20.0, (real_t)1.0)}; // half extents; 2m tall block
        physics::RigidState groundState; groundState.position = Vector3r(0.0, 0.0, (real_t)-1.0); // top face at z=0
        physics::RigidProps groundProps; groundProps.visual = true; groundProps.collidable = true; groundProps.mass = std::nullopt;
        engine.addRigidBody(groundShape, groundState, groundProps);

        const real_t layerWidth = (real_t)0.12; // spacing between leg pair slices
            entt::entity middleAxle = buildWalker(engine, params, L, baseAlphaDeg, nLegPairs, p2_world, planeNormal, xAxis, layerWidth);
        (void)middleAxle; // currently unused in populate
    }

private:
    struct LegGeom2D {
        std::unordered_map<std::string, Vector2r> nodes;
    };

    struct LegGeom3D {
        std::unordered_map<std::string, Vector3r> nodes;
    };

    struct LegBuildResult {
        entt::entity mBeam{entt::null};
        entt::entity foot{entt::null};
        Vector3r p1World{Vector3r::Zero()};
        std::optional<Vector3r> hiWorld;
        std::vector<entt::entity> addedEntities;
    };

    static Vector3r makeOrtho(const Vector3r& n) {
        Vector3r t = (std::abs(n.x()) < 0.5) ? Vector3r::UnitX() : Vector3r::UnitZ();
        Vector3r u = n.cross(t);
        if (u.squaredNorm() < (real_t)1e-10) u = n.cross(Vector3r::UnitY());
        return u.normalized();
    }

    static std::optional<std::pair<Vector2r, Vector2r>> circleIntersections(const Vector2r& c1, real_t r1,
                                                                            const Vector2r& c2, real_t r2) {
        Vector2r d = c2 - c1;
        real_t dist2 = d.squaredNorm();
        real_t dist = std::sqrt((double)dist2);
        if (dist < (real_t)1e-8) return std::nullopt;
        real_t a = (r1*r1 - r2*r2 + dist2) / ((real_t)2 * dist);
        real_t h2 = r1*r1 - a*a;
        if (h2 < (real_t)0) return std::nullopt;
        real_t h = std::sqrt((double)h2);
        Vector2r p = c1 + (a / dist) * d;
        Vector2r perp(-d.y() / dist, d.x() / dist);
        Vector2r p1 = p + h * perp;
        Vector2r p2 = p - h * perp;
        return std::make_pair(p1, p2);
    }

    static std::optional<Vector2r> selectPoint(const std::optional<std::pair<Vector2r, Vector2r>>& pts,
                                               const std::string& mode) {
        if (!pts) return std::nullopt;
        auto [a,b] = *pts;
        if (mode == "upper") return (a.y() >= b.y()) ? a : b;
        if (mode == "lower") return (a.y() <  b.y()) ? a : b;
        if (mode == "left")  return (a.x() <= b.x()) ? a : b;
        if (mode == "right") return (a.x() >  b.x()) ? a : b;
        return a;
    }

    static Vector2r requirePoint(const std::optional<Vector2r>& p, const char* label) {
        if (!p) throw std::runtime_error(std::string("Failed to construct point ") + label);
        return *p;
    }

    // Build leg geometry in 2D (p1 at origin, p2 at (dx,dz))
    static LegGeom2D buildLeg2D(const LegParams& params, const Lengths& L, real_t alpha_deg) {
        const Vector2r p1(0.0, 0.0);
        const Vector2r p2(params.dx, params.dz);
        LegGeom2D leg;

        const real_t alpha = alpha_deg * (real_t)M_PI / (real_t)180.0;
        leg.nodes["mjk"] = p2 + Vector2r(L.m * std::cos(alpha), L.m * std::sin(alpha));
        leg.nodes["ebj"] = requirePoint(selectPoint(circleIntersections(p1, L.b, leg.nodes["mjk"], L.j), "upper"), "ebj");
        leg.nodes["edf"] = requirePoint(selectPoint(circleIntersections(p1, L.d, leg.nodes["ebj"], L.e), "left"), "edf");
        leg.nodes["gcki"] = requirePoint(selectPoint(circleIntersections(p1, L.c, leg.nodes["mjk"], L.k), "lower"), "gcki");
        leg.nodes["fgh"] = requirePoint(selectPoint(circleIntersections(leg.nodes["edf"], L.f, leg.nodes["gcki"], L.g), "left"), "fgh");
        leg.nodes["hi"] = requirePoint(selectPoint(circleIntersections(leg.nodes["fgh"], L.h, leg.nodes["gcki"], L.i), "lower"), "hi");
        leg.nodes["p1"] = p1;
        leg.nodes["p2"] = p2;
        return leg;
    }

    // Build leg geometry in 3D given crank angle and bases.
    LegGeom3D buildLegGeom3D(const LegParams& params,
                             const Lengths& L,
                             real_t alpha_deg,
                             const Vector3r& p2_world,
                             const Vector3r& planeNormal,
                             const Vector3r& xAxis,
                             const std::optional<Vector3r>& mjkOverride3D = std::nullopt) {
        Vector3r n = planeNormal.normalized();
        if (n.squaredNorm() < (real_t)1e-10) n = Vector3r::UnitY();
        Vector3r u = xAxis - xAxis.dot(n) * n;
        if (u.squaredNorm() < (real_t)1e-10) u = makeOrtho(n);
        u.normalize();
        Vector3r v = n.cross(u).normalized();
        // if (v.dot(Vector3r::UnitZ()) < (real_t)0.0) v = -v; // keep vertical orientation consistent across mirrored legs

        Vector3r p1_world = p2_world - (params.dx * u + params.dz * v);

        LegGeom2D leg2d = buildLeg2D(params, L, alpha_deg);
        LegGeom3D leg3d;
        for (const auto& kv : leg2d.nodes) {
            const auto& name = kv.first;
            const Vector2r& pt = kv.second;
            leg3d.nodes[name] = p1_world + pt.x() * u + pt.y() * v;
        }
        if (mjkOverride3D) leg3d.nodes["mjk"] = *mjkOverride3D;
        return leg3d;
    }

    // Build entire leg from precomputed geometry, optionally reusing a shared m beam.
    LegBuildResult buildLeg(physics::PhysicsEngine& engine,
                  LegGeom3D leg3d,
                  entt::entity frameBar,
                  entt::entity sharedM = entt::null,
                  const std::optional<Vector3r>& mjkOverride = std::nullopt,
                  bool createMBeam = true) {
        Vector3r n = Vector3r::UnitY();
        if (mjkOverride) leg3d.nodes["mjk"] = *mjkOverride;

        LegBuildResult out;
        // Create bars
        std::unordered_map<std::string, entt::entity> segEnt;
        const real_t defaultThickness = 0.02;
        const real_t defaultDensity = 800.0;
        std::vector<entt::entity> addedEntities;
        segEnt["frame"] = frameBar;
        if (frameBar != entt::null) addedEntities.push_back(frameBar);
        for (const auto& kv : kSegments) {
            const auto& name = kv.first;
            if (name == "frame") continue; // already created
            const auto& a = kv.second.first;
            const auto& b = kv.second.second;
            if (!leg3d.nodes.count(a) || !leg3d.nodes.count(b)) {
                segEnt[name] = entt::null;
                continue;
            }
            if (name == "m") {
                if (sharedM != entt::null) segEnt[name] = sharedM;
                else if (createMBeam) segEnt[name] = addBar(engine, leg3d.nodes.at(a), leg3d.nodes.at(b), defaultThickness, defaultDensity, /*collidable*/ true);
                else segEnt[name] = entt::null;
            } else {
                segEnt[name] = addBar(engine, leg3d.nodes.at(a), leg3d.nodes.at(b), defaultThickness, defaultDensity, /*collidable*/ true);
            }
            if (segEnt[name] != entt::null) addedEntities.push_back(segEnt[name]);
        }

        // Feet disabled
        entt::entity foot = entt::null;
        out.hiWorld = std::nullopt;

        // Hinge bars at a node using a chain (n-1 hinges instead of all pairs)
        auto hingeAt = [&](const std::string& node, const std::vector<std::string>& bars){
            if (!leg3d.nodes.count(node) || bars.size() < 2) return;
            for (size_t i = 0; i + 1 < bars.size(); ++i) {
                auto ei = segEnt[bars[i]];
                auto ej = segEnt[bars[i + 1]];
                if (ei == entt::null || ej == entt::null) continue;
                physics::JointFrame jf = physics::JointFrame::fromAxis(leg3d.nodes.at(node), n);
                engine.addTranslationRotationConstraint(ei, ej, jf,
                    K_TRANS, D_TRANS, R_ROT, D_ROT);
            }
        };

        hingeAt("p2",    {"frame","m"});
        hingeAt("p1",    {"frame","b","d","c"});
        hingeAt("mjk",   {"m","k","j"});
        hingeAt("ebj",   {"j","b","e"});
        hingeAt("edf",   {"e","f","d"});
        hingeAt("gcki",  {"k","c","g","i"});
        hingeAt("fgh",   {"f","g","h"});
        hingeAt("hi",    {"h","i"});

        last_leg_world_ = leg3d;
        out.mBeam = segEnt.count("m") ? segEnt["m"] : entt::null;
        out.foot = foot;
        if (leg3d.nodes.count("p1")) out.p1World = leg3d.nodes.at("p1");
        out.addedEntities = addedEntities;
        return out;
    }

    // Build two mirrored legs sharing p2, mirrored across plane by flipping xAxis.
    std::pair<LegBuildResult, LegBuildResult> buildLegPair(physics::PhysicsEngine& engine,
                  const LegParams& params,
                  const Lengths& L,
                  real_t alpha_deg,
                  const Vector3r& p2_world,
                  const Vector3r& planeNormal,
                  const Vector3r& xAxis,
                  entt::entity frameBar) {
        LegGeom3D legA = buildLegGeom3D(params, L, alpha_deg, p2_world, -planeNormal, xAxis);
        std::optional<Vector3r> mjkShared = legA.nodes.count("mjk") ? std::make_optional(legA.nodes.at("mjk")) : std::nullopt;
        LegGeom3D legB = buildLegGeom3D(params, L, 180.0- alpha_deg, p2_world, planeNormal, -xAxis);

        // First leg builds the shared m (p2->mjk); second leg reuses the same beam but keeps its own geometry
        entt::entity sharedM = entt::null;

        auto a = buildLeg(engine, legA, frameBar, sharedM, mjkShared, /*createMBeam*/ true);
        sharedM = a.mBeam; // capture created m beam from first leg

        auto b = buildLeg(engine, legB, frameBar, sharedM, std::nullopt, /*createMBeam*/ false);

        a.mBeam = sharedM;
        b.mBeam = sharedM;
        return {a, b};
    }

    // Build walker with phased leg pairs, weld all m-beams together.
    entt::entity buildWalker(physics::PhysicsEngine& engine,
                  const LegParams& params,
                  const Lengths& L,
                  real_t baseAlphaDeg,
                  int nPairs,
                  const Vector3r& p2_world,
                  const Vector3r& planeNormal,
                  const Vector3r& xAxis,
                  real_t layerWidth) {
        std::vector<entt::entity> mBeams;
        std::vector<entt::entity> addedEntities;
        auto pushEntity = [&](entt::entity e){
            if (e == entt::null) return;
            if (std::find(addedEntities.begin(), addedEntities.end(), e) == addedEntities.end()) {
                addedEntities.push_back(e);
            }
        };
        std::vector<Vector3r> p1LeftPoints;
        std::vector<Vector3r> p1RightPoints;
        std::vector<Vector3r> p2Points;
        std::vector<entt::entity> frames;
        std::vector<Vector3r> frameCenters;
        const real_t pairSpacing = layerWidth;
        const real_t frameMargin = (real_t)0.02;
        const real_t axleNormalOffset = (real_t)0.04;
        const real_t walkerLiftLocal = p2_world.z() - params.dz; // derive lift from provided p2
        const real_t axleZMin = params.dz + walkerLiftLocal; // p2 z
        const real_t axleZBase = walkerLiftLocal;            // p1 z
        const real_t axleMin = std::min(axleZBase, axleZMin);
        const real_t axleMax = std::max(axleZBase, axleZMin);
        const real_t halfHeight = (axleMax - axleMin) * (real_t)0.5 + (real_t)0.01; // very tight vertical margin around axles
        const real_t axleCenterZ = (axleMax + axleMin) * (real_t)0.5;
        Vector3r n = planeNormal.normalized();
        Vector3r u = xAxis - xAxis.dot(n) * n; if (u.squaredNorm() < (real_t)1e-10) u = makeOrtho(n); u.normalize();

        // Precompute axle anchor span along u to size frames in x
        real_t minU = std::numeric_limits<real_t>::infinity();
        real_t maxU = -std::numeric_limits<real_t>::infinity();
        for (int i = 0; i < nPairs; ++i) {
            Vector3r p2Offset = p2_world + (real_t)i * pairSpacing * n;
            p2Points.push_back(p2Offset);
            real_t phaseAlpha = baseAlphaDeg + (360.0 / nPairs) * i;
            LegGeom3D legA = buildLegGeom3D(params, L, phaseAlpha, p2Offset, -planeNormal, xAxis);
            LegGeom3D legB = buildLegGeom3D(params, L, 180.0- phaseAlpha, p2Offset, planeNormal, -xAxis);
            if (legA.nodes.count("p1")) {
                real_t proj = (legA.nodes.at("p1") - p2_world).dot(u);
                minU = std::min(minU, proj);
                maxU = std::max(maxU, proj);
            }
            if (legB.nodes.count("p1")) {
                real_t proj = (legB.nodes.at("p1") - p2_world).dot(u);
                minU = std::min(minU, proj);
                maxU = std::max(maxU, proj);
            }
        }
        if (!std::isfinite(minU) || !std::isfinite(maxU)) { minU = -(real_t)0.1; maxU = (real_t)0.1; }
        real_t halfWidthU = std::max(std::abs(minU), std::abs(maxU)) + frameMargin;

        // Frames between pair slices (nPairs-1), rigid and dynamic, oriented with leg plane
        for (int i = 0; i + 1 < nPairs; ++i) {
            Vector3r center = p2_world + ((real_t)i + (real_t)0.5) * pairSpacing * n;
            center.z() = axleCenterZ; // vertically center around axles
            Vector3r v = n.cross(u).normalized();
            Vector3r halfExtents(halfWidthU, // span x to cover axle spacing
                                 (real_t)0.01, // thin along normal between leg planes
                                 halfHeight); // height just to cover axles
            Quaternion4r q; {
                Matrix33r R; R.col(0) = u; R.col(1) = n; R.col(2) = v; q = Quaternion4r(R);
            }
            entt::entity frameBar = addFrameCube(engine, center, halfExtents, q, /*density*/ 500.0);
            frames.push_back(frameBar);
            frameCenters.push_back(center);
            pushEntity(frameBar);
        }


        // Build pairs now that frames exist (each pair uses nearest mid-frame: last pair uses last frame)
        for (int i = 0; i < nPairs; ++i) {
            Vector3r p2Offset = p2_world + (real_t)i * pairSpacing * n;
            size_t fidx = (frames.empty()) ? 0 : std::min<size_t>(i, frames.size() - 1);
            entt::entity frameForPair = frames.empty() ? entt::null : frames[fidx];
            real_t phaseAlpha = baseAlphaDeg + (360.0 / nPairs) * i + 180.0 * i;
            auto pair = buildLegPair(engine, params, L, phaseAlpha, p2Offset, planeNormal, xAxis, frameForPair);
            for (auto e : pair.first.addedEntities) pushEntity(e);
            for (auto e : pair.second.addedEntities) pushEntity(e);
            if (pair.first.mBeam != entt::null) mBeams.push_back(pair.first.mBeam);
            // Capture p1 positions from each leg in the pair
            p1LeftPoints.push_back(pair.first.p1World);
            p1RightPoints.push_back(pair.second.p1World);
        }

        auto addSingleAxle = [&](const std::vector<Vector3r>& pts)->entt::entity{
            if (pts.size() < 2) return entt::null;
            Vector3r start = pts.front();
            Vector3r end   = pts.back();
            return addBar(engine, start, end, (real_t)0.01, (real_t)800.0, /*collidable*/ false);
        };

        entt::entity leftAxle = addSingleAxle(p1LeftPoints); Vector3r leftPos = p1LeftPoints.front() + (p1LeftPoints.back() - p1LeftPoints.front()) * (real_t)0.5;
        entt::entity rightAxle = addSingleAxle(p1RightPoints); Vector3r rightPos = p1RightPoints.front() + (p1RightPoints.back() - p1RightPoints.front()) * (real_t)0.5;
        entt::entity middleAxle = addSingleAxle(p2Points); Vector3r middlePos = p2Points.front() + (p2Points.back() - p2Points.front()) * (real_t)0.5;

        // Attach single axles to all frames via hinge (free about normal, locked otherwise)
        auto attachAxleToFrames = [&](entt::entity axle, Vector3r axlePos = Vector3r::Zero()){
            if (axle == entt::null || frames.empty()) return;
            for (size_t i = 0; i < frames.size(); ++i) {
                if (frames[i] == entt::null) continue;
                physics::JointFrame jf = physics::JointFrame::fromAxis(axlePos, n);
                engine.addTranslationRotationConstraint(axle, frames[i], jf,
                    K_TRANS, D_TRANS, R_ROT, D_ROT);
            }
        };
        attachAxleToFrames(leftAxle, leftPos);
        attachAxleToFrames(rightAxle, rightPos);

        // Track all added entities for global pairwise disabling
        pushEntity(leftAxle);
        pushEntity(rightAxle);
        pushEntity(middleAxle);
        for (auto f : frames) pushEntity(f);
        for (auto m : mBeams) pushEntity(m);

        // Disable collisions between all added entities (axles, frames, bars, m-beams)
        for (size_t i = 0; i < addedEntities.size(); ++i) {
            for (size_t j = i + 1; j < addedEntities.size(); ++j) {
                engine.disableCollisionBetween(addedEntities[i], addedEntities[j]);
            }
        }

        // Attach middle axle rigidly to each m beam at corresponding p2 slice
        if (middleAxle != entt::null) {
            for (size_t i = 0; i < mBeams.size() && i < p2Points.size(); ++i) {
                if (mBeams[i] == entt::null) continue;
                engine.addRigidConstraint(middleAxle, mBeams[i]);
            }
        }

        return middleAxle;
    }

    static entt::entity addBar(physics::PhysicsEngine& engine,
                               const Vector3r& a,
                               const Vector3r& b,
                               real_t thickness,
                               real_t density,
                               bool collidable = true) {
        (void)thickness; (void)density; (void)collidable;
        Vector3r d = b - a;
        real_t len = d.norm();
        if (len < (real_t)1e-6) return entt::null;
        Vector3r mid = (a + b) * (real_t)0.5;
        Quaternion4r q = Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), d.normalized());
        real_t radius = thickness * (real_t)0.5;
        real_t halfLen = len * (real_t)0.5; // extend beyond nodes so the capsule envelopes them
        physics::CapsuleShape shape{radius, halfLen};
        physics::RigidState state; state.position = mid; state.orientation = q;
        physics::RigidProps props;
        props.visual = true;
        props.collidable = collidable;
        if (density > 0) props.density = density; else props.mass = std::nullopt;
        return engine.addRigidBody(shape, state, props);
    }

    static entt::entity addFrameCube(physics::PhysicsEngine& engine,
                                     const Vector3r& center,
                                     const Vector3r& halfExtents,
                                     const Quaternion4r& orientation,
                                     real_t density) {
        (void)density;
        physics::CubeShape shape{halfExtents};
        physics::RigidState state; state.position = center; state.orientation = orientation;
        physics::RigidProps props;
        props.visual = true;
        props.collidable = false;
        if (density > 0) props.density = density; else props.mass = std::nullopt;
        return engine.addRigidBody(shape, state, props);
    }

    std::optional<LegGeom3D> last_leg_world_;
};
