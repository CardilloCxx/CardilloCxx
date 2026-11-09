#pragma once

#include "../net/RingGeometry.hpp"
#include <cardillo/physics/assembly/dynamics_assembler.hpp>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cardillo;

// Pull test on five rockfall net rings, every wire of a ring modelled as a Cosserat rod:
//
//   | O     O |
//   |    O    |  ->
//   | O     O |
//
// Two rings hang on pins of each of two rigid blocks (pins between two plates), and a centre ring links all four diagonally (as
// in a 4-neighbour ring net). The left block is fixed, the right block is pulled in +x (no gravity).
// A ring is a bundle of n closed loops of a steel wire: ROCCO 12/3/300 is a 3 mm wire wound 12 times
// to a ring of 300 mm inner diameter, held together by three clips [1, 2]. The wires interact by
// frictional contact. The pull force on the right pins and the reaction force on the left pins over
// the pin displacement go to <output.folder>/<prefix>_force.csv; cf. the ring tests of Grassl [1] in
// [2], Figs. 10 and 11.
//
// Clips (scene.clip_mode):
//   0  none: the bundle holds together by contact only
//   1  pinned: a rigid clip holds the wire nodes at its position, the wires rotate freely in it
//   2  rigid: the wire nodes are rigidly fixed to the clip (translation and rotation)
//   3  band: a closed steel band (Cosserat rod, rectangular section) is pre-tensioned around the
//      bundle; it clamps the wires by contact pressure and holds them by friction, so they can slip
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Z.-H. Zhu, J.-H. Yin, C.-J. Ouyang, D.-Y. Tan, J.-Q. Qin, Modeling a flexible ring net with the
//     discrete element method, J. Eng. Mech. 146(2), doi:10.1061/(ASCE)EM.1943-7889.0001707.
//
// Scene parameters (scene.<name>, defaults in parentheses):
//   ring_windings (12), wire_diameter (0.003), ring_inner_diameter (0.3)   ROCCO 12/3/300 [1, 2]
//   wire_gap (1e-4)           initial gap between the wires of a bundle (m)
//   wire_segments (48)        beam elements per wire loop
//   link_slack (0.01)         initial slack of the links: centre distance below the taut one (m)
//   crossing_clearance (1.5)  gap of linked rings at their crossings in bundle diameters
//   pin_diameter (0.02), pin_gap (0.002)
//   pull_speed (0.5), pull_ramp (0.02)   speed of the right block (m/s), smooth start over pull_ramp (s)
//   clip_mode (1), clip_count (3), clip_offset_deg (0)   clip angles offset_deg + k 360°/count
//   band_thickness (0.0015), band_width (0.015), band_segments (12), band_pretension (1e-3)
//   primitive_colliders (0)   1: one capsule per wire segment, 0: inter-segment hulls

// Wire of a ring bundle: the bundle centreline offset by (y, z) along its normal and binormal.
class BundleWireSpline : public misc::SplinePattern {
   public:
    BundleWireSpline(const misc::SplinePattern& centreline, const Vector2r& offset) : m_centreline(&centreline), m_offset(offset) {
        constexpr int n = 4096;
        for (int i = 0; i < n; ++i) m_len += (position((real_t)(i + 1) / n) - position((real_t)i / n)).norm();
    }

    real_t totalLength() const override { return m_len; }
    bool isLoop() const override { return true; }
    Vector3r centerOfMass() const override { return Vector3r::Zero(); }

    misc::SplineSample sample(real_t alpha) const override {
        constexpr real_t h = 1e-6;
        const Vector3r tangent = (position(alpha + h) - position(alpha - h)).normalized();
        const misc::SplineSample s = m_centreline->sample(alpha - std::floor(alpha));
        const Vector3r normal = (s.normal - s.normal.dot(tangent) * tangent).normalized();
        return {position(alpha), tangent, normal, tangent.cross(normal)};
    }

   private:
    Vector3r position(real_t alpha) const {
        const misc::SplineSample s = m_centreline->sample(alpha - std::floor(alpha));
        return s.position + m_offset.x() * s.normal + m_offset.y() * s.binormal;
    }

    const misc::SplinePattern* m_centreline{nullptr};
    Vector2r m_offset{Vector2r::Zero()};
    real_t m_len{(real_t)0};
};

class NetDetailedScene : public SceneBase {
   public:
    const char* sceneName() const override { return "net_detailed"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();

        // ---- Wires and bundle: ROCCO 12/3/300 [1, 2]
        const int windings = (int)cfg.sceneParam("ring_windings", 12);
        const real_t wireRadius = cfg.sceneParam("wire_diameter", 0.003) / 2;
        const real_t innerRadius = cfg.sceneParam("ring_inner_diameter", 0.3) / 2;
        const real_t wireGap = cfg.sceneParam("wire_gap", 1e-4);
        const size_t segments = (size_t)cfg.sceneParam("wire_segments", 48);
        const real_t E = 2.0e11;  // high-tensile steel wire [2]
        const real_t nu = 0.3;
        const real_t steelDensity = 7850;

        const std::vector<Vector2r> layout = bundleLayout(windings, 2 * wireRadius + wireGap);
        real_t bundleRadius = 0;
        for (const Vector2r& p : layout) bundleRadius = std::max(bundleRadius, p.norm());
        bundleRadius += wireRadius;
        const real_t ringRadius = innerRadius + bundleRadius;  // centreline of the bundle

        // ---- Rings: centre ring at the origin, corner rings at (±a, ±a), linked with link_slack below
        // the taut centre distance 2 (R - r_b).
        const real_t linkDistance = 2 * (ringRadius - bundleRadius) - cfg.sceneParam("link_slack", 0.01);
        const real_t a = linkDistance / std::sqrt((real_t)2);
        const real_t clearance = cfg.sceneParam("crossing_clearance", 1.5);
        const real_t saddle = saddleAmplitude(ringRadius, Vector3r(a, a, 0), clearance * 2 * bundleRadius);
        const SaddleRingSpline centreline(ringRadius, saddle);
        std::vector<BundleWireSpline> wires;
        for (const Vector2r& p : layout) wires.emplace_back(centreline, p);

        // Initial gaps between the wires of a bundle and between the wires of linked rings.
        constexpr int subdivisions = 8;
        real_t bundleGap = std::numeric_limits<real_t>::infinity();
        real_t linkGap = std::numeric_limits<real_t>::infinity();
        for (size_t i = 0; i < wires.size(); ++i) {
            for (size_t j = i + 1; j < wires.size(); ++j) bundleGap = std::min(bundleGap, polygonDistance(wires[i], Vector3r::Zero(), wires[j], Vector3r::Zero(), segments, subdivisions));
            for (size_t j = 0; j < wires.size(); ++j) linkGap = std::min(linkGap, polygonDistance(wires[i], Vector3r::Zero(), wires[j], Vector3r(a, a, 0), segments, subdivisions));
        }
        bundleGap -= 2 * wireRadius;
        linkGap -= 2 * wireRadius;
        if (bundleGap <= 0 || linkGap <= 0) throw std::runtime_error("NetDetailedScene: wires intersect initially; increase scene.wire_gap, scene.crossing_clearance or scene.wire_segments");

        const bool primitive = cfg.sceneParam("primitive_colliders", 0) != 0;
        const auto colliderMode = primitive ? physics::BeamColliderMode::RigidBodyPrimitive : physics::BeamColliderMode::InterSegmentHull;
        const auto wireSection = physics::BeamCrossSection::round(wireRadius);
        const auto wireSprings = physics::BeamSpringParams::fromMaterial(E, nu);
        const auto wireProps = physics::RigidProps::withDensity(steelDensity);

        // ---- Clips
        const int clipMode = (int)cfg.sceneParam("clip_mode", 1);
        const int clipCount = (int)cfg.sceneParam("clip_count", 3);
        const size_t clipOffset = (size_t)std::lround(cfg.sceneParam("clip_offset_deg", 0) / 360 * (real_t)segments) % segments;
        const real_t bandThickness = cfg.sceneParam("band_thickness", 0.0015);
        const real_t bandWidth = cfg.sceneParam("band_width", 0.015);
        const size_t bandSegments = (size_t)cfg.sceneParam("band_segments", 12);
        const real_t bandRadius = bundleRadius + wireGap + bandThickness / 2;
        const real_t clipMass = steelDensity * (real_t)(2 * M_PI) * bandRadius * bandThickness * bandWidth;
        const auto bandSection = physics::BeamCrossSection::rectangular(bandThickness, bandWidth);
        auto bandSprings = physics::BeamSpringParams::fromMaterial(E, nu);
        bandSprings.withAxialPretension(cfg.sceneParam("band_pretension", 1e-3));
        // Beam nodes sit at alpha = k / segments (hulls) or at the segment midpoints (capsules).
        const real_t nodeShift = primitive ? (real_t)0.5 : (real_t)0;

        const std::array<Vector3r, 5> centres{Vector3r::Zero(), Vector3r(-a, a, 0), Vector3r(-a, -a, 0), Vector3r(a, a, 0), Vector3r(a, -a, 0)};
        for (const Vector3r& centre : centres) {
            std::vector<std::vector<entt::entity>> nodes;  // per wire, in the order of alpha
            for (const BundleWireSpline& wire : wires) {
                const entt::entity root = engine.createBeam(wire, wireSection, wireSprings, physics::RigidState(centre, Quaternion4r::Identity()), wireProps, segments, colliderMode).first;
                nodes.push_back(chain(engine, root));
            }
            if (clipMode == 0) continue;
            for (int c = 0; c < clipCount; ++c) {
                const size_t k = (clipOffset + (size_t)c * segments / (size_t)clipCount) % segments;
                const misc::SplineSample s = centreline.sample(((real_t)k + nodeShift) / (real_t)segments);
                const Vector3r position = centre + s.position;
                if (clipMode == 3) {
                    const misc::CircleSpline band(position, bandRadius, s.tangent, s.normal);
                    engine.createBeam(band, bandSection, bandSprings, physics::RigidState(Vector3r::Zero(), Quaternion4r::Identity()), physics::RigidProps::withDensity(steelDensity), bandSegments,
                                      colliderMode);
                    continue;
                }
                Matrix33r R;
                R << s.tangent, s.normal, s.binormal;
                const entt::entity clip = engine.addRigidBody(physics::CubeShape(Vector3r(bandWidth / 2, bandRadius, bandRadius)), physics::RigidState(position, Quaternion4r(R)),
                                                              physics::RigidProps(clipMass, -1, true, false));
                for (const auto& wireNodes : nodes) {
                    if (clipMode == 1) {
                        engine.addTranslationalConstraint(clip, wireNodes[k], physics::JointFrame(wireNodes[k]));
                    } else {
                        engine.addRigidConstraint(clip, wireNodes[k]);
                    }
                }
            }
        }

        // ---- Pins along z through the outer side of the corner rings, each held between two plates
        // above and below the ring (the blocks). Left: fixed, right: pulled in +x.
        const real_t pinRadius = cfg.sceneParam("pin_diameter", 0.02) / 2;
        const real_t pinGap = cfg.sceneParam("pin_gap", 0.002);
        const real_t pinX = a + ringRadius - bundleRadius - pinRadius - pinGap;
        m_pullSpeed = cfg.sceneParam("pull_speed", 0.5);
        m_pullRamp = cfg.sceneParam("pull_ramp", 0.02);
        // The plates cover the bundle around the pin; their faces lie pin_gap above and below the
        // bundle within that footprint (corner ring at (-a, a), pin at (-pinX, a); all corners alike).
        const real_t plateHalfWidth = 2 * (pinRadius + pinGap + bundleRadius);
        real_t zLow = std::numeric_limits<real_t>::infinity();
        real_t zHigh = -std::numeric_limits<real_t>::infinity();
        constexpr int samples = 4096;
        for (int i = 0; i < samples; ++i) {
            const Vector3r q = Vector3r(-a, a, 0) + centreline.sample((real_t)i / samples).position;
            if (std::abs(q.x() + pinX) > plateHalfWidth + bundleRadius || std::abs(q.y() - a) > plateHalfWidth + bundleRadius) continue;
            zLow = std::min(zLow, q.z() - bundleRadius - pinGap);
            zHigh = std::max(zHigh, q.z() + bundleRadius + pinGap);
        }
        const real_t pinZ = (zLow + zHigh) / 2;
        const real_t pinHalfLength = (zHigh - zLow) / 2;
        const Vector3r plateHalfExtents(plateHalfWidth, plateHalfWidth, pinRadius);
        for (const real_t side : {(real_t)-1, (real_t)1}) {
            std::vector<std::pair<entt::entity, Vector3r>> bodies;
            for (const real_t y : {a, -a}) {
                const Vector3r position(side * pinX, y, pinZ);
                const entt::entity pin = engine.addRigidBody(physics::CapsuleShape(pinRadius, pinHalfLength), physics::RigidState(position), physics::RigidProps(0));
                (side > 0 ? m_rightPins : m_leftPins).push_back(pin);
                bodies.emplace_back(pin, position);
                for (const real_t z : {zLow - plateHalfExtents.z(), zHigh + plateHalfExtents.z()}) {
                    const Vector3r platePosition(side * pinX, y, z);
                    bodies.emplace_back(engine.addRigidBody(physics::CubeShape(plateHalfExtents), physics::RigidState(platePosition), physics::RigidProps(0)), platePosition);
                }
            }
            if (side < 0) continue;
            for (const auto& [body, position] : bodies) {
                const real_t v = m_pullSpeed;
                const real_t ramp = m_pullRamp;
                engine.addTrajectory(body, [position, v, ramp](real_t t) -> TrajectoryPose { return TrajectoryPose{position + Vector3r::UnitX() * displacement(v, ramp, t), Quaternion4r::Identity()}; });
            }
        }

        std::filesystem::create_directories(cfg.output_folder);
        m_forceFile.open(cfg.output_folder + "/" + cfg.outputPrefix() + "_force.csv");
        m_forceFile << "t,displacement,pull_force,left_force\n";

        const real_t wireMass = steelDensity * (real_t)M_PI * wireRadius * wireRadius * centreline.totalLength();
        std::cout << "Rings: 5 x " << windings << " wires (" << windings * wireMass << " kg per ring), bundle radius " << bundleRadius << " m, gaps: bundle " << bundleGap << " m, link "
                  << linkGap << " m\n"
                  << "  link distance " << linkDistance << " m, saddle amplitude " << saddle << " m, clip mode " << clipMode << "\n";
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t) override {
        // Contact forces of the last step on the pins (the impulse acts on body b of a contact).
        Vector3r right = Vector3r::Zero();
        Vector3r left = Vector3r::Zero();
        for (const auto& c : engine.dynamicsAssembler().contacts()) {
            const Vector3r impulse = c.last_impulse.x() * c.normal + c.last_impulse.y() * c.tangent1 + c.last_impulse.z() * c.tangent2;
            auto add = [&](const std::vector<entt::entity>& pins, Vector3r& sum) {
                if (std::find(pins.begin(), pins.end(), c.b) != pins.end()) sum += impulse;
                if (std::find(pins.begin(), pins.end(), c.a) != pins.end()) sum -= impulse;
            };
            add(m_rightPins, right);
            add(m_leftPins, left);
        }
        const real_t h = engine.timeStep();
        m_forceFile << t << ',' << displacement(m_pullSpeed, m_pullRamp, t) << ',' << -right.x() / h << ',' << left.x() / h << '\n';
    }

   private:
    // Displacement of the right block: speed v reached smoothly (constant acceleration) over the ramp time.
    static real_t displacement(real_t v, real_t ramp, real_t t) { return t < ramp ? v * t * t / (2 * ramp) : v * (t - ramp / 2); }

    // Cross-section positions of n wires with centre distance `spacing`: the n points of a hexagonal
    // lattice closest to a lattice point or to a lattice triangle centre, whichever is more compact.
    static std::vector<Vector2r> bundleLayout(int n, real_t spacing) {
        const Vector2r e1(spacing, 0);
        const Vector2r e2(spacing / 2, spacing * std::sqrt((real_t)3) / 2);
        std::vector<Vector2r> best;
        real_t bestRadius = std::numeric_limits<real_t>::infinity();
        for (const Vector2r& origin : std::array<Vector2r, 2>{Vector2r::Zero(), (e1 + e2) / 3}) {
            std::vector<Vector2r> points;
            for (int i = -n; i <= n; ++i)
                for (int j = -n; j <= n; ++j) points.push_back(i * e1 + j * e2 - origin);
            std::sort(points.begin(), points.end(), [](const Vector2r& p, const Vector2r& q) { return p.squaredNorm() < q.squaredNorm(); });
            points.resize((size_t)n);
            if (points.back().norm() < bestRadius) {
                bestRadius = points.back().norm();
                best = points;
            }
        }
        return best;
    }

    // Nodes of the closed beam starting at `root`, in order.
    static std::vector<entt::entity> chain(physics::PhysicsEngine& engine, entt::entity root) {
        std::vector<entt::entity> nodes;
        entt::entity e = root;
        do {
            nodes.push_back(e);
            e = engine.ecs().get<C_BeamElement>(e).next.value_or(root);
        } while (e != root);
        return nodes;
    }

    std::vector<entt::entity> m_rightPins;
    std::vector<entt::entity> m_leftPins;
    real_t m_pullSpeed{(real_t)0};
    real_t m_pullRamp{(real_t)0};
    std::ofstream m_forceFile;
};
