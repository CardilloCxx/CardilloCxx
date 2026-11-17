#pragma once

#include "RingGeometry.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

// Rings modelled wire by wire (NetDetailedScene, NetWiresScene). ROCCO n/d/D denotes a ring of a d mm
// steel wire wound n times to D mm inner diameter, held together by clips [1, 2]. Every wire is a
// Cosserat rod; the wires interact by frictional contact. Ring models (scene.ring_wire):
//   0  n closed loops, one per winding
//   1  one wire wound n times through n + 1 slots of the bundle cross-section, each pass changing to the
//      next slot in its own short window of the turn (woundChain(), WoundWireSpline); both ends extend
//      wire_overhang past the first clip, so that they overlap and the clip holds them
//
// Clips (scene.clip_mode), clip_count per ring at the beam nodes nearest to the angles
// clip_offset_deg + k 360° / clip_count:
//   0  none: the bundle holds together by contact only
//   1  pinned: a rigid steel sleeve (cylinder around the bundle) holds the wire nodes at its position,
//      the wires rotate freely in it
//   2  rigid: the wire nodes are rigidly fixed to the sleeve (translation and rotation)
//   3  band: a closed steel band (Cosserat rod, rectangular section) is pre-tensioned around the
//      bundle; it clamps the wires by contact pressure and holds them by friction, so they can slip
// Sleeves (band_width long, band_thickness wall) and bands collide with the other rings, the shackles
// and the block; a sleeve does not collide with the wire elements it encloses.
//
// Scene parameters (scene.<name>, defaults in parentheses; the scenes set the defaults marked *):
//   ring_wire (*), ring_windings (*), wire_diameter (0.003), ring_inner_diameter (0.3)
//   wire_gap (1e-4)           initial gap between the wires (passes) of a bundle (m)
//   wire_segments (*)         beam elements per turn
//   wire_overhang (0.03)      length by which each end of the wound wire extends past the centre of the
//                             first clip (m), rounded up to whole beam elements of the flat ring
//   clip_mode (1), clip_count (3), clip_offset_deg (0)
//   band_thickness (0.0015), band_width (0.015)   steel of the sleeves and bands (m)
//   band_segments (12), band_pretension (1e-3)     band elements, relative shortening of the band
//   primitive_colliders (0)   1: one capsule per wire segment, 0: inter-segment hulls
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Z.-H. Zhu, J.-H. Yin, C.-J. Ouyang, D.-Y. Tan, J.-Q. Qin, Modeling a flexible ring net with the
//     discrete element method, J. Eng. Mech. 146(2), doi:10.1061/(ASCE)EM.1943-7889.0001707.

// Cross-section positions of n wires with centre distance `spacing`: the n points of a hexagonal
// lattice closest to a lattice point or to a lattice triangle centre, whichever is more compact. A
// partly filled outermost shell is filled evenly over its angles.
inline std::vector<Vector2r> bundleLayout(int n, real_t spacing) {
    const Vector2r e1(spacing, 0);
    const Vector2r e2(spacing / 2, spacing * std::sqrt((real_t)3) / 2);
    // Lattice points relative to the origin as (A e1 + B e2) / 3 with integers A = 3i - o, B = 3j - o
    // (o = 0: lattice point, o = 1: triangle centre); 9 |p|² / spacing² = A² + AB + B² is exact.
    struct Point {
        long A{0};
        long B{0};
        long key{0};
    };
    std::vector<Point> best;
    for (const long o : {0L, 1L}) {
        std::vector<Point> points;
        for (long i = -n; i <= n; ++i) {
            for (long j = -n; j <= n; ++j) {
                const long A = 3 * i - o;
                const long B = 3 * j - o;
                points.push_back({A, B, A * A + A * B + B * B});
            }
        }
        std::sort(points.begin(), points.end(), [](const Point& p, const Point& q) { return p.key < q.key; });
        points.resize((size_t)n);
        if (best.empty() || points.back().key < best.back().key) best = points;
    }
    auto position = [&](const Point& p) -> Vector2r { return ((real_t)p.A * e1 + (real_t)p.B * e2) / 3; };

    // Outermost shell: all lattice points at the distance of the last chosen one, chosen evenly by angle.
    const long o = (best.front().A % 3 + 3) % 3 == 0 ? 0 : 1;
    const long outerKey = best.back().key;
    std::vector<Vector2r> shell;
    for (long i = -n; i <= n; ++i) {
        for (long j = -n; j <= n; ++j) {
            const Point p{3 * i - o, 3 * j - o, (3 * i - o) * (3 * i - o) + (3 * i - o) * (3 * j - o) + (3 * j - o) * (3 * j - o)};
            if (p.key == outerKey) shell.push_back(position(p));
        }
    }
    std::sort(shell.begin(), shell.end(), [](const Vector2r& p, const Vector2r& q) { return std::atan2(p.y(), p.x()) < std::atan2(q.y(), q.x()); });
    std::vector<Vector2r> layout;
    for (const Point& p : best)
        if (p.key < outerKey) layout.push_back(position(p));
    const size_t needed = (size_t)n - layout.size();
    for (size_t k = 0; k < needed; ++k) layout.push_back(shell[k * shell.size() / needed]);
    return layout;
}

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

// Slots of a wire wound n times (scene.ring_wire = 1) in the bundle cross-section: n + 1 points of
// bundleLayout(n + 1) in chain order. Pass k lies in slot k and changes to slot k + 1 in its own short
// window of the turn (WoundWireSpline), while all other passes rest in their slots; the moving pass
// comes closest to the slots other than k and k + 1. The chain is the order that keeps that distance
// largest (depth-first search over all orders, pruned by the best order found so far; for a compact
// cluster a chain of neighbouring slots, sqrt(3)/2 of their spacing), scaled so that the passes are at
// least `spacing` apart.
inline std::vector<Vector2r> woundChain(int n, real_t spacing) {
    const std::vector<Vector2r> points = bundleLayout(n + 1, 1);
    std::vector<size_t> best;
    real_t bestDistance = 0;
    std::vector<size_t> chain;
    std::vector<bool> used(points.size(), false);
    // Closest approach of the move from slot a to slot b to the slots not yet visited and visited before
    // a, i.e. all slots other than a and b.
    auto moveDistance = [&](size_t a, size_t b) {
        real_t d = std::numeric_limits<real_t>::infinity();
        const Vector2r step = points[b] - points[a];
        for (size_t i = 0; i < points.size(); ++i) {
            if (i == a || i == b) continue;
            const real_t f = std::clamp((points[i] - points[a]).dot(step) / step.squaredNorm(), (real_t)0, (real_t)1);
            d = std::min(d, (points[a] + f * step - points[i]).norm());
        }
        return d;
    };
    auto extend = [&](auto&& self, real_t distance) -> void {
        if (chain.size() == points.size()) {
            best = chain;
            bestDistance = distance;
            return;
        }
        for (size_t i = 0; i < points.size(); ++i) {
            if (used[i]) continue;
            const real_t d = chain.empty() ? distance : std::min(distance, moveDistance(chain.back(), i));
            if (d <= bestDistance) continue;
            used[i] = true;
            chain.push_back(i);
            self(self, d);
            chain.pop_back();
            used[i] = false;
        }
    };
    extend(extend, std::numeric_limits<real_t>::infinity());
    std::vector<Vector2r> slots;
    for (const size_t i : best) slots.push_back(points[i] * spacing / bestDistance);
    return slots;
}

// Wire wound n times around the bundle centreline through the slots of woundChain(), from the ring
// parameter alpha0 - overhang to alpha0 + n + overhang (u in [-overhang, n + overhang] turns from alpha0):
// both ends overlap by `overhang` turns past alpha0 (the first clip), the start in slot 0, the end in
// slot n. In turn u = k + f (f in [0, 1)) pass k lies in slot k and moves smoothly to slot k + 1 in its
// own window; the n windows tile f in [z, 1 - z] in reverse order (the last pass moves first, into the
// free slot n, every other pass into the slot just vacated). All passes rest for |f| < z around alpha0,
// z = max(1 / (2 (n + 1)), overhang), so the overlapping ends lie beside resting passes only.
class WoundWireSpline : public misc::SplinePattern {
   public:
    WoundWireSpline(const misc::SplinePattern& centreline, const std::vector<Vector2r>& slots, real_t alpha0, real_t overhang)
        : m_centreline(&centreline), m_slots(slots), m_alpha0(alpha0), m_overhang(overhang) {
        const int n = turns();
        m_rest = std::max((real_t)0.5 / (real_t)(n + 1), overhang);
        m_window = (1 - 2 * m_rest) / (real_t)n;
        constexpr int samplesPerTurn = 4096;
        const int samples = (int)std::ceil(samplesPerTurn * span());
        for (int i = 0; i < samples; ++i) m_len += (position(-m_overhang + span() * (real_t)(i + 1) / samples) - position(-m_overhang + span() * (real_t)i / samples)).norm();
    }

    int turns() const { return (int)m_slots.size() - 1; }
    // Turns from one end to the other.
    real_t span() const { return (real_t)turns() + 2 * m_overhang; }
    real_t totalLength() const override { return m_len; }
    bool isLoop() const override { return false; }
    Vector3r centerOfMass() const override { return Vector3r::Zero(); }

    misc::SplineSample sample(real_t alpha) const override {
        constexpr real_t h = 1e-6;
        const real_t u = -m_overhang + span() * std::clamp(alpha, (real_t)0, (real_t)1);
        const real_t u0 = std::max(u - h, -m_overhang);
        const real_t u1 = std::min(u + h, (real_t)turns() + m_overhang);
        const Vector3r tangent = (position(u1) - position(u0)).normalized();
        const misc::SplineSample s = m_centreline->sample(ringAlpha(u));
        const Vector3r normal = (s.normal - s.normal.dot(tangent) * tangent).normalized();
        return {position(u), tangent, normal, tangent.cross(normal)};
    }

   private:
    real_t ringAlpha(real_t u) const { return m_alpha0 + u - std::floor(m_alpha0 + u); }

    Vector3r position(real_t u) const {
        const int n = turns();
        const int k = std::clamp((int)std::floor(u), 0, n - 1);
        const real_t x = std::clamp((u - k - m_rest - (real_t)(n - 1 - k) * m_window) / m_window, (real_t)0, (real_t)1);
        const Vector2r offset = m_slots[k] + x * x * (3 - 2 * x) * (m_slots[k + 1] - m_slots[k]);
        const misc::SplineSample s = m_centreline->sample(ringAlpha(u));
        return s.position + offset.x() * s.normal + offset.y() * s.binormal;
    }

    const misc::SplinePattern* m_centreline{nullptr};
    std::vector<Vector2r> m_slots;
    real_t m_alpha0{(real_t)0};
    real_t m_overhang{(real_t)0};
    real_t m_rest{(real_t)0};    // half-width of the zone around alpha0 where all passes rest (turns)
    real_t m_window{(real_t)0};  // width of a slot-change window (turns)
    real_t m_len{(real_t)0};
};

// Wires, material and clips of a wire-resolved ring, read from the scene parameters.
struct WireRingModel {
    static constexpr real_t E = 2.0e11;  // high-tensile steel wire [2]
    static constexpr real_t nu = 0.3;
    static constexpr real_t steelDensity = 7850;

    using Wires = std::vector<std::unique_ptr<misc::SplinePattern>>;
    using RingNodes = std::vector<std::vector<entt::entity>>;  // beam nodes per wire, in the order of the wire

    bool wound{false};
    int windings{0};
    real_t wireRadius{0};
    real_t wireGap{0};
    size_t segments{0};            // beam elements per turn
    size_t endElements{0};         // beam elements by which each end of the wound wire extends past the first clip
    real_t bundleRadius{0};        // envelope of the wires
    real_t ringRadius{0};          // centreline of the bundle
    std::vector<Vector2r> layout;  // cross-section positions of the loops, or the slots of the wound wire
    physics::BeamColliderMode colliderMode{physics::BeamColliderMode::InterSegmentHull};
    real_t nodeShift{0};  // beam nodes sit at the element ends (hulls) or at their midpoints (capsules)

    int clipMode{0};
    int clipCount{0};
    size_t clipOffset{0};  // node index of the first clip
    real_t bandThickness{0};
    real_t bandWidth{0};
    size_t bandSegments{0};
    real_t bandRadius{0};  // centreline of the band
    real_t bandPretension{0};

    WireRingModel(const config::Config& cfg, int defaultWire, int defaultWindings, size_t defaultSegments) {
        wound = cfg.sceneParam("ring_wire", (real_t)defaultWire) != 0;
        windings = (int)cfg.sceneParam("ring_windings", (real_t)defaultWindings);
        wireRadius = cfg.sceneParam("wire_diameter", 0.003) / 2;
        wireGap = cfg.sceneParam("wire_gap", 1e-4);
        segments = (size_t)cfg.sceneParam("wire_segments", (real_t)defaultSegments);
        layout = wound ? woundChain(windings, 2 * wireRadius + wireGap) : bundleLayout(windings, 2 * wireRadius + wireGap);
        for (const Vector2r& p : layout) bundleRadius = std::max(bundleRadius, p.norm());
        bundleRadius += wireRadius;
        ringRadius = cfg.sceneParam("ring_inner_diameter", 0.3) / 2 + bundleRadius;
        if (wound) {
            const real_t elementLength = (real_t)(2 * M_PI) * ringRadius / (real_t)segments;
            endElements = (size_t)std::ceil(cfg.sceneParam("wire_overhang", 0.03) / elementLength);
            if (2 * endElements >= segments) throw std::runtime_error("WireRingModel: scene.wire_overhang must be shorter than half a turn");
        }
        const bool primitive = cfg.sceneParam("primitive_colliders", 0) != 0;
        colliderMode = primitive ? physics::BeamColliderMode::RigidBodyPrimitive : physics::BeamColliderMode::InterSegmentHull;
        nodeShift = primitive ? (real_t)0.5 : (real_t)0;

        clipMode = (int)cfg.sceneParam("clip_mode", 1);
        clipCount = (int)cfg.sceneParam("clip_count", 3);
        clipOffset = (size_t)std::lround(cfg.sceneParam("clip_offset_deg", 0) / 360 * (real_t)segments) % segments;
        bandThickness = cfg.sceneParam("band_thickness", 0.0015);
        bandWidth = cfg.sceneParam("band_width", 0.015);
        bandSegments = (size_t)cfg.sceneParam("band_segments", 12);
        bandRadius = bundleRadius + wireGap + bandThickness / 2;
        bandPretension = cfg.sceneParam("band_pretension", 1e-3);
    }

    // Beam elements of one wire.
    size_t wireSegments() const { return wound ? (size_t)windings * segments + 2 * endElements : segments; }
    // Ring parameter of ring node k (k = 0..segments-1).
    real_t nodeAlpha(size_t k) const { return ((real_t)k + nodeShift) / (real_t)segments; }
    // Ring node index of clip c.
    size_t clipNode(int c) const { return (clipOffset + (size_t)c * segments / (size_t)clipCount) % segments; }
    // Ring node index nearest to the ring parameter alpha.
    size_t nearestNode(real_t alpha) const { return (size_t)((std::lround(alpha * (real_t)segments - nodeShift) % (long)segments + (long)segments) % (long)segments); }
    // Outer radius of a clip (sleeve or band).
    real_t clipRadius() const { return bandRadius + bandThickness / 2; }

    // Wires of a ring with the given centreline (ring centre at the origin). The ends of the wound wire
    // overlap at the first clip.
    Wires wires(const misc::SplinePattern& centreline) const {
        Wires result;
        if (wound) {
            result.push_back(std::make_unique<WoundWireSpline>(centreline, layout, (real_t)clipOffset / (real_t)segments, (real_t)endElements / (real_t)segments));
        } else {
            for (const Vector2r& p : layout) result.push_back(std::make_unique<BundleWireSpline>(centreline, p));
        }
        return result;
    }

    // Beam nodes of all wires (passes) of a ring at ring node index k.
    std::vector<entt::entity> nodesAt(const RingNodes& nodes, size_t k) const {
        std::vector<entt::entity> result;
        if (!wound) {
            for (const auto& wireNodes : nodes) result.push_back(wireNodes[k]);
            return result;
        }
        // Node i of the wound wire lies at ring node clipOffset - endElements + i.
        const size_t rel = (k + segments - clipOffset + endElements) % segments;
        for (size_t i = rel; i < nodes.front().size(); i += segments) result.push_back(nodes.front()[i]);
        return result;
    }

    // Smallest initial gap between the wire surfaces within a bundle: between different loops, or between
    // passes of the wound wire at least half a turn apart along it.
    real_t bundleGap(const Wires& wires) const {
        real_t d = std::numeric_limits<real_t>::infinity();
        if (wound) {
            const std::vector<Vector3r> points = polygon(*wires.front(), Vector3r::Zero());
            const size_t halfTurn = segments * (size_t)subdivisions / 2;
            for (size_t i = 0; i < points.size(); ++i)
                for (size_t j = i + halfTurn; j < points.size(); ++j) d = std::min(d, (points[j] - points[i]).norm());
        } else {
            for (size_t i = 0; i < wires.size(); ++i)
                for (size_t j = i + 1; j < wires.size(); ++j) d = std::min(d, polygonDistance(*wires[i], Vector3r::Zero(), *wires[j], Vector3r::Zero(), segments, subdivisions));
        }
        return d - 2 * wireRadius;
    }

    // Smallest initial gap between the wire surfaces of two rings `offset` apart.
    real_t linkGap(const Wires& wires, const Vector3r& offset) const {
        real_t d = std::numeric_limits<real_t>::infinity();
        for (const auto& a : wires)
            for (const auto& b : wires) d = std::min(d, polygonDistance(*a, Vector3r::Zero(), *b, offset, wireSegments(), subdivisions));
        return d - 2 * wireRadius;
    }

    // Smallest initial gap between the clips of a ring and the wires of rings at `offsets`: distance of
    // the wire centrelines from the cylinder of radius clipRadius() and length band_width around the
    // bundle, minus the wire radius.
    real_t clipGap(const misc::SplinePattern& centreline, const Wires& wires, const std::vector<Vector3r>& offsets) const {
        real_t gap = std::numeric_limits<real_t>::infinity();
        if (clipMode == 0) return gap;
        for (int c = 0; c < clipCount; ++c) {
            const misc::SplineSample s = centreline.sample(nodeAlpha(clipNode(c)));
            for (const Vector3r& offset : offsets) {
                for (const auto& wire : wires) {
                    for (const Vector3r& point : polygon(*wire, offset)) {
                        const Vector3r x = point - s.position;
                        const real_t axial = x.dot(s.tangent);
                        const real_t radial = (x - axial * s.tangent).norm();
                        const Vector2r outside(std::max(std::abs(axial) - bandWidth / 2, (real_t)0), std::max(radial - clipRadius(), (real_t)0));
                        gap = std::min(gap, outside.norm() - wireRadius);
                    }
                }
            }
        }
        return gap;
    }

    // Wires and clips of a ring at `centre`; returns the beam nodes of every wire.
    RingNodes addRing(physics::PhysicsEngine& engine, const misc::SplinePattern& centreline, const Wires& wires, const Vector3r& centre) const {
        const auto wireSection = physics::BeamCrossSection::round(wireRadius);
        const auto wireSprings = physics::BeamSpringParams::fromMaterial(E, nu);
        const auto wireProps = physics::RigidProps::withDensity(steelDensity);
        RingNodes nodes;
        for (const auto& wire : wires) {
            const entt::entity root = engine.createBeam(*wire, wireSection, wireSprings, physics::RigidState(centre, Quaternion4r::Identity()), wireProps, wireSegments(), colliderMode).first;
            nodes.push_back(beamNodes(engine, root));
        }
        if (clipMode == 0) return nodes;

        const auto bandSection = physics::BeamCrossSection::rectangular(bandThickness, bandWidth);
        auto bandSprings = physics::BeamSpringParams::fromMaterial(E, nu);
        bandSprings.withAxialPretension(bandPretension);
        const real_t clipMass = steelDensity * (real_t)(2 * M_PI) * bandRadius * bandThickness * bandWidth;
        // Wire elements a sleeve encloses: those within half its width (plus the wire radius) of the
        // clipped node, counted with a lower bound of the element length (inner wire of a flat ring).
        const real_t elementLength = (real_t)(2 * M_PI) * (ringRadius - bundleRadius) / (real_t)segments;
        const size_t sleeveReach = (size_t)std::ceil((bandWidth / 2 + wireRadius) / elementLength);
        for (int c = 0; c < clipCount; ++c) {
            const size_t k = clipNode(c);
            const misc::SplineSample s = centreline.sample(nodeAlpha(k));
            const Vector3r position = centre + s.position;
            if (clipMode == 3) {
                const misc::CircleSpline band(position, bandRadius, s.tangent, s.normal);
                engine.createBeam(band, bandSection, bandSprings, physics::RigidState(Vector3r::Zero(), Quaternion4r::Identity()), physics::RigidProps::withDensity(steelDensity), bandSegments,
                                  colliderMode);
                continue;
            }
            // Sleeve: steel tube around the bundle, a cylinder with its axis (local z) along the bundle.
            Matrix33r R;
            R << s.normal, s.binormal, s.tangent;
            const entt::entity clip = engine.addRigidBody(physics::CylinderShape(clipRadius(), bandWidth / 2), physics::RigidState(position, Quaternion4r(R)), physics::RigidProps(clipMass));
            for (const entt::entity node : nodesAt(nodes, k)) {
                if (clipMode == 1) {
                    engine.addTranslationalConstraint(clip, node, physics::JointFrame(node));
                } else {
                    engine.addRigidConstraint(clip, node);
                }
                disableAlongWire(engine, clip, node, sleeveReach);
            }
        }
        return nodes;
    }

    // Steel mass of the wires of a ring.
    real_t ringMass(const Wires& wires) const {
        real_t length = 0;
        for (const auto& wire : wires) length += wire->totalLength();
        return steelDensity * (real_t)M_PI * wireRadius * wireRadius * length;
    }

   private:
    static constexpr int subdivisions = 8;  // samples per beam element in the clearance checks

    // Disable the contacts of `body` with the wire node `node` and the `reach` nodes before and after it
    // along the wire. A node carries the collider of its element (hull to the next node, or capsule).
    static void disableAlongWire(physics::PhysicsEngine& engine, entt::entity body, entt::entity node, size_t reach) {
        engine.disableCollisionBetween(body, node);
        for (const bool forward : {true, false}) {
            entt::entity e = node;
            for (size_t i = 0; i < reach; ++i) {
                const auto& element = engine.ecs().get<C_BeamElement>(e);
                const std::optional<entt::entity> neighbour = forward ? element.next : element.prev;
                if (!neighbour) break;
                e = *neighbour;
                engine.disableCollisionBetween(body, e);
            }
        }
    }

    // Points along the polygon through the beam nodes of a wire, shifted by `offset`.
    std::vector<Vector3r> polygon(const misc::SplinePattern& wire, const Vector3r& offset) const {
        const size_t n = wireSegments();
        std::vector<Vector3r> points;
        for (size_t k = 0; k < n; ++k) {
            const Vector3r p0 = wire.sample((real_t)k / (real_t)n).position + offset;
            const Vector3r p1 = wire.sample((real_t)(k + 1) / (real_t)n).position + offset;
            for (int s = 0; s < subdivisions; ++s) points.push_back(p0 + (p1 - p0) * ((real_t)s / subdivisions));
        }
        return points;
    }
};
