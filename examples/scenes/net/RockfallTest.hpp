#pragma once

#include "RingGeometry.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cardillo;

// Test rig of the rockfall ring net scenes (NetScene, NetWiresScene): the frame the net is shackled to,
// the shackles, and the falling block with optional debris (gravity -z, net in the x-y plane at z = 0).
//
// Scene parameters (scene.<name>, defaults in parentheses):
//   block_mass (850), block_reference_density (2700), block_mesh (0: rock.obj, 1: rock2.obj)
//   drop_height (3), start_gap (= drop_height)   fall height and initial gap above the net in m
//   extra_boulders (0), extra_size_min (0.15), extra_size_max (0.4), extra_density (2700),
//   extra_spread (2), extra_seed (1)    additional smaller rocks falling with the block
//
// References
// [1] Eberhardt et al., 1D elastoplastic models with isotropic hardening as differential inclusions in
//     primal form, JTCAM, 2026, Sec. 4.2 (supplementary/ring_nets/Eberhardt2026.pdf).
// [2] EOTA, ETAG 027: Guideline for European technical approval of falling rock protection kits,
//     April 2013, Annex A.2, Figure A.2.

// Rectangular frame around the net, visual only; the shackle strings end on its inner edges, the anchor
// lines x = ±anchorX and y = ±anchorY in the plane z = 0.
class NetFrame {
   public:
    NetFrame(physics::PhysicsEngine& engine, real_t anchorX, real_t anchorY) {
        const real_t profile = 0.2;
        for (size_t i = 0; i < m_sides.size(); ++i) {
            const Vector3r& d = m_sides[i];
            const bool alongY = d.x() != 0;
            const real_t anchorLine = alongY ? anchorX : anchorY;
            const real_t halfSpan = (alongY ? anchorY : anchorX) + 2 * profile;
            const Vector3r halfExtents = alongY ? Vector3r(profile / 2, halfSpan, profile / 2) : Vector3r(halfSpan, profile / 2, profile / 2);
            m_positions[i] = d * (anchorLine + profile / 2) - Vector3r::UnitZ() * (profile / 2);
            m_bodies[i] = engine.addRigidBody(physics::CubeShape(halfExtents), physics::RigidState(m_positions[i]), physics::RigidProps(0, -1, true, false));
        }
    }

    // Elastic string from the point rLocal (body coordinates) of `body` to `anchor` on the frame side
    // facing the outward direction `d`; rest length = initial length.
    void addString(physics::PhysicsEngine& engine, entt::entity body, const Vector3r& d, const Vector3r& anchor, real_t stiffness, real_t damping, const Vector3r& rLocal = Vector3r::Zero()) const {
        const size_t i = (size_t)(std::find(m_sides.begin(), m_sides.end(), d) - m_sides.begin());
        if (i == m_sides.size()) throw std::runtime_error("NetFrame: shackle direction must be ±x or ±y");
        engine.addLinearDistanceConstraint(m_bodies[i], body, anchor - m_positions[i], rLocal, stiffness, damping);
    }

   private:
    std::array<Vector3r, 4> m_sides{Vector3r::UnitX(), -Vector3r::UnitX(), Vector3r::UnitY(), -Vector3r::UnitY()};
    std::array<entt::entity, 4> m_bodies{entt::null, entt::null, entt::null, entt::null};
    std::array<Vector3r, 4> m_positions{Vector3r::Zero(), Vector3r::Zero(), Vector3r::Zero(), Vector3r::Zero()};
};

// Shackle of a boundary ring (scene.shackle_mode = 1): a rigid steel ring (torus) around the ring
// bundle, hung on the frame by an elastic string from its outermost point. The shackle is far stiffer
// than the bundle; it holds the bundle by contact only, so the ring can slide in it. The torus mesh is
// written to <output.folder>/<prefix>_shackle.obj.
//
// Scene parameters (scene.<name>, defaults in parentheses):
//   shackle_diameter (0.012)   diameter of the shackle's steel section (m)
//   shackle_gap (0.003)        initial gap between the shackle and the bundle envelope (m)
//   shackle_resolution (24)    facets of the torus mesh around the ring (half as many around the section)
class ShackleModel {
   public:
    ShackleModel(const config::Config& cfg, real_t bundleRadius) {
        m_rodRadius = cfg.sceneParam("shackle_diameter", 0.012) / 2;
        const real_t gap = cfg.sceneParam("shackle_gap", 0.003);
        const int resolution = (int)cfg.sceneParam("shackle_resolution", 24);
        // The facets lie inside the torus, so the mesh clears the bundle by at least `gap`.
        m_radius = bundleRadius + gap + m_rodRadius;
        std::filesystem::create_directories(cfg.output_folder);
        m_mesh = cfg.output_folder + "/" + cfg.outputPrefix() + "_shackle.obj";
        writeTorus(m_mesh, m_radius, m_rodRadius, resolution, std::max(resolution / 2, 3));
    }

    real_t rodRadius() const { return m_rodRadius; }
    real_t radius() const { return m_radius; }  // centreline radius of the torus

    // Shackle around the bundle at `p` (bundle tangent `t`), hung by a string of rest length
    // stringLength - radius() to the frame anchor p + stringLength d; returns the shackle.
    entt::entity add(physics::PhysicsEngine& engine, const NetFrame& frame, const Vector3r& p, const Vector3r& t, const Vector3r& d, real_t stringLength, real_t stiffness, real_t damping) const {
        if (stringLength <= m_radius) throw std::runtime_error("ShackleModel: scene.string_length must exceed the shackle radius");
        const entt::entity shackle =
            engine.addRigidBody(physics::MeshShape(m_mesh), physics::RigidState(p, Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), t)), physics::RigidProps::withDensity(steelDensity));
        // Attachment: the outermost point of the torus centreline, in the body frame (principal axes).
        const VectorXr state = engine.getPosition(shackle);
        Quaternion4r q;
        q.coeffs() = state.tail<4>();
        const Vector3r rLocal = q.toRotationMatrix().transpose() * (p + m_radius * d - state.head<3>());
        frame.addString(engine, shackle, d, p + stringLength * d, stiffness, damping, rLocal);
        return shackle;
    }

   private:
    static constexpr real_t steelDensity = 7850;

    // Closed torus about the z-axis: centreline radius R, section radius r, nu x nv quads as triangles.
    static void writeTorus(const std::string& path, real_t R, real_t r, int nu, int nv) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("ShackleModel: cannot write " + path);
        f << "# torus R = " << R << " m, r = " << r << " m (NetScene / NetWiresScene shackle)\n";
        for (int i = 0; i < nu; ++i) {
            const real_t u = (real_t)(2 * M_PI) * i / nu;
            for (int j = 0; j < nv; ++j) {
                const real_t v = (real_t)(2 * M_PI) * j / nv;
                f << "v " << (R + r * std::cos(v)) * std::cos(u) << ' ' << (R + r * std::cos(v)) * std::sin(u) << ' ' << r * std::sin(v) << '\n';
            }
        }
        auto index = [&](int i, int j) { return ((i + nu) % nu) * nv + (j + nv) % nv + 1; };
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                f << "f " << index(i, j) << ' ' << index(i + 1, j) << ' ' << index(i + 1, j + 1) << '\n';
                f << "f " << index(i, j) << ' ' << index(i + 1, j + 1) << ' ' << index(i, j + 1) << '\n';
            }
        }
    }

    real_t m_rodRadius{0};
    real_t m_radius{0};
    std::string m_mesh;
};

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
inline MeshMoments meshMoments(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("meshMoments: cannot open " + path);
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

// Falling block above a net whose top lies at z = netTop, and optional debris; returns the block.
inline entt::entity addRockfallBlock(physics::PhysicsEngine& engine, real_t netTop) {
    const auto& cfg = engine.config();
    const std::string meshDir = std::string(CARDILLO_SOURCE_DIR) + "/res/meshes/";

    // Block: a non-convex rock with the mass and inertia of the EOTA block [2] of [1] (850 kg).
    // The EOTA block is a cube of edge L with square faces of edge L/2 and cut edges and corners, of
    // density 2500-3000 kg/m³ [2]; its volume is 17/24 L³ and its second moments per unit volume are
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

    // [1]: dropped from h = 3 m, impact velocity 7.67 m/s, no rotation. The fall height is measured
    // from the block's bounding sphere to the top of the net; start_gap < drop_height starts the
    // block closer to the net with the corresponding velocity (the net then starts unsagged).
    const real_t dropHeight = cfg.sceneParam("drop_height", 3);
    const real_t startGap = std::min(cfg.sceneParam("start_gap", dropHeight), dropHeight);
    const real_t g = engine.gravity().norm();
    const Vector3r v0 = -Vector3r::UnitZ() * std::sqrt(2 * g * (dropHeight - startGap));
    const Vector3r blockPosition = Vector3r::UnitZ() * (netTop + startGap + blockRadius);
    physics::RigidProps blockProps;
    blockProps.mass = blockMass;
    const entt::entity block = engine.addRigidBody(physics::MeshShape(rockMeshes[blockMesh], blockScale), physics::RigidState(blockPosition, v0, Quaternion4r::Identity()), blockProps);
    engine.track(block, "block");

    // Optional debris: smaller non-convex rocks of random size (bounding diameter) and orientation in
    // layers above the block, falling with it.
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
        const Vector3r position((col - (perRow - 1) * (real_t)0.5) * sizeMax, (row - (perRow - 1) * (real_t)0.5) * sizeMax, blockPosition.z() + blockRadius + (layer + (real_t)0.5) * sizeMax);
        const size_t mesh = (size_t)k % rockMeshes.size();
        const real_t size = sizeMin + (sizeMax - sizeMin) * uniform(rng);
        const real_t scale = size / (2 * rockMoments[mesh].boundingRadius(Vector3r::Ones()));
        const Quaternion4r q = Quaternion4r(normal(rng), normal(rng), normal(rng), normal(rng)).normalized();
        const entt::entity rock =
            engine.addRigidBody(physics::MeshShape(rockMeshes[mesh], Vector3r::Constant(scale)), physics::RigidState(position, v0, q), physics::RigidProps::withDensity(extraDensity));
        extraMass += engine.getMass(rock)(0, 0);
    }

    std::cout << "Block: " << blockMass << " kg, EOTA reference edge " << L << " m, mesh scale " << blockScale.transpose() << ", density " << blockMass / (bm.volume * blockScale.prod())
              << " kg/m^3, inertia " << engine.getMass(block).diagonal().tail<3>().transpose() << " kg m^2 (EOTA " << blockMass * (real_t)42 / 340 * L * L << ")\n"
              << "  impact velocity " << std::sqrt(2 * g * dropHeight) << " m/s, energy " << blockMass * g * dropHeight << " J\n";
    if (extraCount > 0) std::cout << "Debris: " << extraCount << " rocks, " << extraMass << " kg\n";
    return block;
}
