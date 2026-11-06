#pragma once

#include "RingGeometry.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cardillo;

// Drop test on a horizontal rockfall ring net (Geobrugg ROCCO rings). This is the field test of
// Grassl [1] (setup in [3], Fig. 12), which Eberhardt et al. [2] use as benchmark (Sec. 4.2,
// Fig. 19c): a square ring net of 3.9 m x 3.9 m in which every ring is linked with four neighbours is
// shackled to the top of a rigid steel frame, and a block of about 850 kg falls 3 m onto it (gravity -z).
// The shackles are modelled as elastic strings (LinearDistanceConstraint) to a static frame.
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Eberhardt et al., 1D elastoplastic models with isotropic hardening as differential inclusions in
//     primal form, JTCAM, 2026, Sec. 4.2 (supplementary/ring_nets/Eberhardt2026.pdf).
// [3] Z.-H. Zhu, J.-H. Yin, C.-J. Ouyang, D.-Y. Tan, J.-Q. Qin, Modeling a flexible ring net with the
//     discrete element method, J. Eng. Mech. 146(2), doi:10.1061/(ASCE)EM.1943-7889.0001707.
// [4] EOTA, ETAG 027: Guideline for European technical approval of falling rock protection kits,
//     April 2013, Annex A.2, Figure A.2.
// [5] Geobrugg, product data sheet RXE-8000 (BAFU): ROCCO ring net 19/3/300, wire diameter 3 mm, wire
//     tensile strength 1770 N/mm². ROCCO n/d/D denotes rings of n windings of a d mm wire, diameter D mm.
//
// Scene parameters (scene.<name>, defaults in parentheses):
//   rings_x, rings_y (10, 10)       rows of rings_x rings alternate with rows of rings_x - 1 rings;
//                                   rings_y long rows, i.e. 10 x 10 + 9 x 9 = 181 rings
//   net_size_x, net_size_y (3.9)    outer net dimensions in m [1, 3]
//   ring_windings (7), wire_diameter (0.003), ring_inner_diameter (0.3)   ROCCO 7/3/300 [2, 3, 5]
//   ring_segments (16)              beam elements per ring
//   axial_stiffness_factor (0.2)    measured / theoretical axial ring stiffness [1, 3]
//   crossing_clearance (1.5)        gap of linked rings at their crossings in bundle diameters
//   string_length (0.1), string_stiffness (1e7), string_damping (0)   shackle strings (m, N/m, Ns/m)
//   block_mass (850), block_reference_density (2700), block_mesh (0: rock.obj, 1: rock2.obj)
//   drop_height (3), start_gap (= drop_height)   fall height and initial gap above the net in m
//   extra_boulders (0), extra_size_min (0.15), extra_size_max (0.4), extra_density (2700),
//   extra_spread (2), extra_seed (1)    additional smaller rocks falling with the block
//   primitive_colliders (0)         1: one capsule per ring segment, 0: inter-segment hulls

class NetScene : public SceneBase {
   public:
    const char* sceneName() const override { return "net"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();
        const std::string meshDir = std::string(CARDILLO_SOURCE_DIR) + "/res/meshes/";

        // ---- Rings: ROCCO 7/3/300, 7 windings of a 3 mm wire [2], 300 mm inner diameter, held by three clips [3]
        const real_t windings = cfg.sceneParam("ring_windings", 7);
        const real_t wireRadius = cfg.sceneParam("wire_diameter", 0.003) / 2;
        const real_t innerRadius = cfg.sceneParam("ring_inner_diameter", 0.3) / 2;
        const size_t segments = (size_t)cfg.sceneParam("ring_segments", 16);

        // The windings form a loose bundle. Contact uses its envelope (a hexagonally packed 7-wire bundle
        // is 3 wire diameters wide, fill factor 7/9); mass and stiffness follow the wires: area n A_w and,
        // for bending and torsion, n I_w of separate wires instead of the solid envelope [3, Eq. 10].
        const real_t bundleRadius = 3 * wireRadius * std::sqrt(windings / 7);
        const real_t ringRadius = innerRadius + bundleRadius;  // centreline
        const real_t areaRatio = windings * wireRadius * wireRadius / (bundleRadius * bundleRadius);
        const real_t inertiaRatio = windings * std::pow(wireRadius / bundleRadius, 4);
        // Grassl measured an axial ring stiffness of 1/5 of n E A_w for rings with 7, 12 and 19 windings,
        // as the windings are not loaded evenly [1], [3, Eq. 14].
        const real_t axialFactor = cfg.sceneParam("axial_stiffness_factor", 0.2);
        const real_t E = 2.0e11;  // high-tensile steel wire [3]; tensile strength 1770 N/mm² [5] (no plasticity here)
        const real_t nu = 0.3;
        const real_t steelDensity = 7850;
        const auto section = physics::BeamCrossSection::round(bundleRadius);
        const auto springs = physics::BeamSpringParams::fromMaterial(E, nu, axialFactor * areaRatio, areaRatio, inertiaRatio, inertiaRatio, inertiaRatio);
        const auto ringProps = physics::RigidProps::withDensity(steelDensity * areaRatio);

        // ---- Net [1, 3]: rows of nx and nx - 1 rings alternate, every ring is linked with its four
        // diagonal neighbours at the offsets (±ax, ±ay).
        const int nx = (int)cfg.sceneParam("rings_x", 10);
        const int ny = (int)cfg.sceneParam("rings_y", 10);
        const real_t netSizeX = cfg.sceneParam("net_size_x", 3.9);
        const real_t netSizeY = cfg.sceneParam("net_size_y", 3.9);
        const real_t outerRadius = ringRadius + bundleRadius;
        const real_t ax = (netSizeX - 2 * outerRadius) / (2 * (nx - 1));
        const real_t ay = (netSizeY - 2 * outerRadius) / (2 * (ny - 1));
        const real_t linkDistance = std::hypot(ax, ay);
        // Linked rings are taut at a centre distance 2 (R - r_b); rings of one row (or column) must not touch.
        if (linkDistance >= 2 * (ringRadius - bundleRadius) || std::min(ax, ay) <= outerRadius) {
            throw std::runtime_error("NetScene: net size does not fit the number and diameter of the rings");
        }

        // Saddle rings (RingGeometry.hpp) start the woven net flat and interlinked; at the crossings the
        // centrelines of linked rings are crossing_clearance bundle diameters apart.
        const real_t clearance = cfg.sceneParam("crossing_clearance", 1.5);
        const real_t saddle = saddleAmplitude(ringRadius, Vector3r(ax, ay, 0), clearance * 2 * bundleRadius);
        const SaddleRingSpline ring(ringRadius, saddle);
        const real_t linkGap = polygonDistance(ring, Vector3r::Zero(), ring, Vector3r(ax, ay, 0), segments) - 2 * bundleRadius;
        if (linkGap <= 0) throw std::runtime_error("NetScene: linked rings intersect initially; increase scene.crossing_clearance or scene.ring_segments");

        const auto colliderMode = cfg.sceneParam("primitive_colliders", 0) != 0 ? physics::BeamColliderMode::RigidBodyPrimitive : physics::BeamColliderMode::InterSegmentHull;

        struct Attachment {
            entt::entity ring{entt::null};
            Vector3r direction{Vector3r::Zero()};
        };
        std::vector<Attachment> attachments;
        int ringCount = 0;
        for (int j = 0; j < 2 * ny - 1; ++j) {
            const bool longRow = j % 2 == 0;
            const int count = longRow ? nx : nx - 1;
            for (int i = 0; i < count; ++i) {
                const Vector3r center((real_t)(2 * i - (count - 1)) * ax, (real_t)(j - (ny - 1)) * ay, 0);
                const entt::entity root = engine.createBeam(ring, section, springs, physics::RigidState(center, Quaternion4r::Identity()), ringProps, segments, colliderMode).first;
                ++ringCount;
                // Shackled boundary rings: both ends of the long rows and the first and last long row.
                if (!longRow) continue;
                if (i == 0) attachments.push_back({root, -Vector3r::UnitX()});
                if (i == count - 1) attachments.push_back({root, Vector3r::UnitX()});
                if (j == 0) attachments.push_back({root, -Vector3r::UnitY()});
                if (j == 2 * ny - 2) attachments.push_back({root, Vector3r::UnitY()});
            }
        }

        // ---- Frame and shackles: each boundary ring hangs on the frame by an elastic string from its
        // outermost node, corner rings on both sides. The frame profile is visual only.
        const real_t stringLength = cfg.sceneParam("string_length", 0.1);
        const real_t stringStiffness = cfg.sceneParam("string_stiffness", 1e7);
        const real_t stringDamping = cfg.sceneParam("string_damping", 0);
        const real_t profile = 0.2;
        const real_t anchorX = (nx - 1) * ax + ringRadius + stringLength;
        const real_t anchorY = (ny - 1) * ay + ringRadius + stringLength;
        auto addFrameSide = [&](const Vector3r& d) {
            const bool alongY = d.x() != 0;
            const real_t anchorLine = alongY ? anchorX : anchorY;
            const real_t halfSpan = (alongY ? anchorY : anchorX) + 2 * profile;
            const Vector3r halfExtents = alongY ? Vector3r(profile / 2, halfSpan, profile / 2) : Vector3r(halfSpan, profile / 2, profile / 2);
            const Vector3r position = d * (anchorLine + profile / 2) - Vector3r::UnitZ() * (profile / 2);
            return std::make_pair(engine.addRigidBody(physics::CubeShape(halfExtents), physics::RigidState(position), physics::RigidProps(0, -1, true, false)), position);
        };
        const std::vector<Vector3r> sides{Vector3r::UnitX(), -Vector3r::UnitX(), Vector3r::UnitY(), -Vector3r::UnitY()};
        std::vector<std::pair<entt::entity, Vector3r>> frame;
        for (const Vector3r& d : sides) frame.push_back(addFrameSide(d));
        for (const Attachment& a : attachments) {
            const size_t side = (size_t)(std::find(sides.begin(), sides.end(), a.direction) - sides.begin());
            const entt::entity node = outermostNode(engine, a.ring, a.direction);
            const Vector3r anchor = engine.getPosition(node).head<3>() + stringLength * a.direction;
            engine.addLinearDistanceConstraint(frame[side].first, node, anchor - frame[side].second, Vector3r::Zero(), stringStiffness, stringDamping);
        }

        // ---- Block: a non-convex rock with the mass and inertia of the EOTA block [4] of [2] (850 kg).
        // The EOTA block is a cube of edge L with square faces of edge L/2 and cut edges and corners, of
        // density 2500-3000 kg/m³ [4]; its volume is 17/24 L³ and its second moments per unit volume are
        // isotropic, <x²> = 21/340 L². The rock mesh is scaled per axis to these second moments.
        const real_t blockMass = cfg.sceneParam("block_mass", 850);
        const real_t referenceDensity = cfg.sceneParam("block_reference_density", 2700);
        const real_t L = std::cbrt(blockMass / referenceDensity * 24 / 17);
        const std::vector<std::string> rockMeshes{meshDir + "rock.obj", meshDir + "rock2.obj"};
        std::vector<MeshMoments> rockMoments;
        for (const auto& path : rockMeshes) rockMoments.push_back(meshMoments(path));
        const size_t blockMesh = cfg.sceneParam("block_mesh", 0) != 0 ? 1 : 0;
        const MeshMoments& bm = rockMoments[blockMesh];
        const Vector3r blockScale = (Vector3r::Constant((real_t)21 / 340 * L * L).cwiseQuotient(bm.secondMoments)).cwiseSqrt();
        // The principal axes of the (nearly isotropic) block are not those of the mesh file, so the
        // bounding sphere about the centroid defines the fall height.
        const real_t blockRadius = bm.boundingRadius(blockScale);

        // [2]: dropped from h = 3 m, impact velocity 7.67 m/s, no rotation. The fall height is measured
        // from the block's bounding sphere to the top of the net; start_gap < drop_height starts the
        // block closer to the net with the corresponding velocity (the net then starts unsagged).
        const real_t dropHeight = cfg.sceneParam("drop_height", 3);
        const real_t startGap = std::min(cfg.sceneParam("start_gap", dropHeight), dropHeight);
        const real_t g = engine.gravity().norm();
        const Vector3r v0 = -Vector3r::UnitZ() * std::sqrt(2 * g * (dropHeight - startGap));
        const real_t netTop = saddle + bundleRadius;
        const Vector3r blockPosition = Vector3r::UnitZ() * (netTop + startGap + blockRadius);
        physics::RigidProps blockProps;
        blockProps.mass = blockMass;
        const entt::entity block = engine.addRigidBody(physics::MeshShape(rockMeshes[blockMesh], blockScale), physics::RigidState(blockPosition, v0, Quaternion4r::Identity()), blockProps);
        engine.track(block, "block");

        // ---- Optional debris: smaller non-convex rocks of random size (bounding diameter) and orientation
        // in layers above the block, falling with it.
        const int extraCount = (int)cfg.sceneParam("extra_boulders", 0);
        const real_t sizeMin = cfg.sceneParam("extra_size_min", 0.15);
        const real_t sizeMax = std::max(cfg.sceneParam("extra_size_max", 0.4), sizeMin);
        const real_t extraDensity = cfg.sceneParam("extra_density", 2700);
        const real_t spread = cfg.sceneParam("extra_spread", 2);
        std::mt19937 rng((unsigned)cfg.sceneParam("extra_seed", 1));
        std::uniform_real_distribution<real_t> uniform(0, 1);
        std::normal_distribution<real_t> normal(0, 1);
        const int perRow = std::max(1, (int)std::floor(spread / sizeMax));
        real_t extraMass = 0;
        for (int k = 0; k < extraCount; ++k) {
            const int col = k % perRow;
            const int row = (k / perRow) % perRow;
            const int layer = k / (perRow * perRow);
            const Vector3r position((col - (perRow - 1) * (real_t)0.5) * sizeMax, (row - (perRow - 1) * (real_t)0.5) * sizeMax,
                                    blockPosition.z() + blockRadius + (layer + (real_t)0.5) * sizeMax);
            const size_t mesh = (size_t)k % rockMeshes.size();
            const real_t size = sizeMin + (sizeMax - sizeMin) * uniform(rng);
            const real_t scale = size / (2 * rockMoments[mesh].boundingRadius(Vector3r::Ones()));
            const Quaternion4r q = Quaternion4r(normal(rng), normal(rng), normal(rng), normal(rng)).normalized();
            const entt::entity rock = engine.addRigidBody(physics::MeshShape(rockMeshes[mesh], Vector3r::Constant(scale)), physics::RigidState(position, v0, q),
                                                          physics::RigidProps::withDensity(extraDensity));
            extraMass += engine.getMass(rock)(0, 0);
        }

        const real_t ringMass = steelDensity * windings * (real_t)M_PI * wireRadius * wireRadius * ring.totalLength();
        std::cout << "Ring net: " << ringCount << " rings (mass " << ringMass << " kg each, " << ringCount * ringMass << " kg), " << attachments.size() << " shackle strings\n"
                  << "  link distance " << linkDistance << " m (taut at " << 2 * (ringRadius - bundleRadius) << " m), saddle amplitude " << saddle << " m, link gap "
                  << linkGap << " m\n"
                  << "Block: " << blockMass << " kg, EOTA reference edge " << L << " m, mesh scale " << blockScale.transpose() << ", density "
                  << blockMass / (bm.volume * blockScale.prod()) << " kg/m^3, inertia " << engine.getMass(block).diagonal().tail<3>().transpose() << " kg m^2 (EOTA "
                  << blockMass * (real_t)42 / 340 * L * L << ")\n"
                  << "  impact velocity " << std::sqrt(2 * g * dropHeight) << " m/s, energy " << blockMass * g * dropHeight << " J\n";
        if (extraCount > 0) std::cout << "Debris: " << extraCount << " rocks, " << extraMass << " kg\n";
    }

   private:
    struct MeshMoments {
        real_t volume{0};
        Vector3r secondMoments{Vector3r::Zero()};  // <x²>, <y²>, <z²> about the centroid per unit volume
        std::vector<Vector3r> vertices;            // relative to the centroid

        // Largest vertex distance from the centroid of the mesh scaled by `scale`.
        real_t boundingRadius(const Vector3r& scale) const {
            real_t r = 0;
            for (const Vector3r& v : vertices) r = std::max(r, scale.cwiseProduct(v).norm());
            return r;
        }
    };

    // Volume moments of a closed triangulated OBJ mesh (fan triangulation of polygons).
    static MeshMoments meshMoments(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("NetScene: cannot open " + path);
        std::vector<Vector3r> vertices;
        real_t volume = 0;
        Vector3r first = Vector3r::Zero();
        Matrix33r second = Matrix33r::Zero();
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream ss(line);
            std::string tag;
            ss >> tag;
            if (tag == "v") {
                Vector3r v;
                ss >> v.x() >> v.y() >> v.z();
                vertices.push_back(v);
            } else if (tag == "f") {
                std::vector<int> face;
                for (std::string token; ss >> token;) face.push_back(std::stoi(token.substr(0, token.find('/'))) - 1);
                for (size_t k = 1; k + 1 < face.size(); ++k) {
                    const Vector3r& a = vertices[face[0]];
                    const Vector3r& b = vertices[face[k]];
                    const Vector3r& c = vertices[face[k + 1]];
                    const real_t v = a.dot(b.cross(c)) / 6;  // tetrahedron with the origin
                    const Vector3r s = a + b + c;
                    volume += v;
                    first += v * s / 4;
                    second += v / 20 * (a * a.transpose() + b * b.transpose() + c * c.transpose() + s * s.transpose());
                }
            }
        }
        MeshMoments m;
        m.volume = volume;
        const Vector3r centroid = first / volume;
        m.secondMoments = (second / volume - centroid * centroid.transpose()).diagonal();
        for (const Vector3r& v : vertices) m.vertices.push_back(v - centroid);
        return m;
    }

    // Node of the ring starting at `root` that lies farthest in direction `d`.
    static entt::entity outermostNode(physics::PhysicsEngine& engine, entt::entity root, const Vector3r& d) {
        entt::entity best = root;
        real_t bestDistance = -std::numeric_limits<real_t>::infinity();
        entt::entity e = root;
        do {
            const real_t distance = engine.getPosition(e).head<3>().dot(d);
            if (distance > bestDistance) {
                bestDistance = distance;
                best = e;
            }
            e = engine.ecs().get<C_BeamElement>(e).next.value_or(root);
        } while (e != root);
        return best;
    }
};
