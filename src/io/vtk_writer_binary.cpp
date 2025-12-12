#include "vtk_writer_binary.hpp"
#include "../collision/collision_coal.hpp"
#include "vtk_sphere_util.hpp"
#include "../solver/warmstart.hpp"
#include "../physics/constraints.hpp"
#include <cmath>
#include <coal/hfield.h>
#include <mpi.h>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include "../partitioning/naive_partitioner.hpp"

namespace fs = std::filesystem;

namespace cardillo::io {

VtkWriterBinary::VtkWriterBinary(std::string outputDir, std::string baseName, int frequency)
    : m_outputDir(std::move(outputDir)), m_baseName(std::move(baseName)), m_frequency(frequency) {
    if (m_frequency < 1) m_frequency = 1;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
}

void VtkWriterBinary::setOutputDir(const std::string& dir) {
    m_outputDir = dir; if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
}
void VtkWriterBinary::setBaseName(const std::string& name) { m_baseName = name; }
void VtkWriterBinary::setFrequency(int freq) { m_frequency = std::max(1, freq); }

void VtkWriterBinary::maybeWrite(int step, real_t time, const cardillo::PhysicsSystem& sys) {
    if (step % m_frequency == 0) write(step, time, sys);
}

void VtkWriterBinary::write(int step, real_t /*time*/, const cardillo::PhysicsSystem& sys) {
    auto sc = sys.timings().scope(cardillo::misc::TimingManager::TimerId::OutputWrite);
    Collected data = collect(sys);
    writePointsOnly(step, 0, data);
    // Write static geometry once or when structure changed
    if (!m_staticGeoWritten || sys.isStructureDirty()) {
        writeStaticGeometry(data);
        m_staticGeoWritten = true;
    }
    // Additionally write a timestep-qualified static geometry snapshot every time
    // writeStaticGeometryStep(step, data);
    // Always write dynamic geometry per step
    writeDynamicGeometry(step, data);
    if (m_writeContacts) {
        try {
            auto& mgr = const_cast<cardillo::PhysicsSystem&>(sys).collisionManager();
            if (sys.consumeStructureDirty()) mgr.rebuild();
            const auto& contacts = mgr.lastFlattenedContacts();
            const bool writeBody = sys.config().output_contacts_body_vectors;
            writeContacts(step, contacts, writeBody, sys.warmstartProvider());
        } catch (...) {
            // skip if collision manager not available
        }
    }
    if (m_writeSprings) {
        writeSprings(step, sys);
    }
}

static inline float f32(real_t v) { return static_cast<float>(v); }

// Partition helper (same heuristic as ASCII writer)
static inline float partitionFromBodyIndex_(const cardillo::PhysicsSystem& sys, const entt::registry& reg, entt::entity e) {
    if (reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(e)) {
        int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(e).b;
        int worldRank = 0, worldSize = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        int Nb = sys.numBodies();
        if (Nb <= 0 || worldSize <= 0) {
            return -1.f; // no valid partition when no bodies or invalid world size
        }
        int bodiesPerRank = std::max(1, (Nb + worldSize - 1) / worldSize);
        return (float)cardillo::partitioning::NaivePartitioner::ownerOf(b, worldSize, bodiesPerRank);
    }
    return -1.f;
}

namespace {

void generateCapsuleMesh(int segmentsCircumference,
                         int segmentsHemisphere,
                         int segmentsCylinder,
                         real_t radius,
                         real_t halfLength,
                         std::vector<Vector3r>& vertices,
                         std::vector<Eigen::Vector3i>& triangles)
{
    vertices.clear();
    triangles.clear();
    if (radius <= (real_t)0) return;
    const int nCirc = std::max(segmentsCircumference, 3);
    const int nHem = std::max(segmentsHemisphere, 1);
    const int nCyl = std::max(segmentsCylinder, 1);

    const real_t pi = std::acos(static_cast<real_t>(-1));
    const real_t halfPi = pi * (real_t)0.5;

    auto addRing = [&](real_t z, real_t ringRadius) {
        std::vector<int> idx;
        idx.reserve(nCirc);
        for (int i = 0; i < nCirc; ++i) {
            const real_t angle = static_cast<real_t>(2.0 * pi * i / nCirc);
            const real_t cx = ringRadius * std::cos(angle);
            const real_t cy = ringRadius * std::sin(angle);
            idx.push_back(static_cast<int>(vertices.size()));
            vertices.emplace_back(cx, cy, z);
        }
        return idx;
    };

    std::vector<std::vector<int>> rings;
    rings.reserve(nHem * 2 + nCyl + 2);

    const int bottomPole = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, -halfLength - radius);

    for (int i = 1; i < nHem; ++i) {
        const real_t theta = (static_cast<real_t>(i) / nHem) * halfPi;
        const real_t ringRadius = radius * std::sin(theta);
        const real_t z = -halfLength - radius * std::cos(theta);
        rings.push_back(addRing(z, ringRadius));
    }

    rings.push_back(addRing(-halfLength, radius));

    for (int i = 1; i < nCyl; ++i) {
        const real_t t = static_cast<real_t>(i) / nCyl;
        const real_t z = -halfLength + t * (halfLength * 2);
        rings.push_back(addRing(z, radius));
    }

    rings.push_back(addRing(halfLength, radius));

    for (int i = nHem - 1; i >= 1; --i) {
        const real_t theta = (static_cast<real_t>(i) / nHem) * halfPi;
        const real_t ringRadius = radius * std::sin(theta);
        const real_t z = halfLength + radius * std::cos(theta);
        rings.push_back(addRing(z, ringRadius));
    }

    const int topPole = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, halfLength + radius);

    if (rings.empty()) return;

    const auto& firstRing = rings.front();
    for (int j = 0; j < nCirc; ++j) {
        const int a = bottomPole;
        const int b = firstRing[j];
        const int c = firstRing[(j + 1) % nCirc];
        triangles.emplace_back(a, b, c);
    }

    for (std::size_t i = 0; i + 1 < rings.size(); ++i) {
        const auto& r0 = rings[i];
        const auto& r1 = rings[i + 1];
        for (int j = 0; j < nCirc; ++j) {
            const int jn = (j + 1) % nCirc;
            triangles.emplace_back(r0[j], r1[j], r0[jn]);
            triangles.emplace_back(r0[jn], r1[j], r1[jn]);
        }
    }

    const auto& lastRing = rings.back();
    for (int j = 0; j < nCirc; ++j) {
        const int a = lastRing[j];
        const int b = lastRing[(j + 1) % nCirc];
        const int c = topPole;
        triangles.emplace_back(a, b, c);
    }
}

} // namespace

VtkWriterBinary::Collected VtkWriterBinary::collect(const cardillo::PhysicsSystem& sys) const {
    Collected out;
    const auto& reg = sys.ecs();

    // Points
    {
        auto vpoints = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                                cardillo::PhysicsSystem::C_PointVisualTag,
                                cardillo::PhysicsSystem::C_Position3,
                                cardillo::PhysicsSystem::C_Mass,
                                cardillo::PhysicsSystem::C_LinearVelocity3,
                                cardillo::PhysicsSystem::C_Radius>();
        for (auto [e, pos, m, v, r] : vpoints.each()) {
            // Skip rigid-body spheres: they will be emitted as triangle meshes
            if (reg.any_of<cardillo::PhysicsSystem::C_RB_Sphere>(e)) continue;
            PointOut po;
            po.pos = pos.value;
            po.mass = static_cast<float>(m.m);
            po.vel = v.value;
            po.radius = static_cast<float>(r.r);
            po.partition = partitionFromBodyIndex_(sys, reg, e);
            po.entityId = static_cast<int>(entt::to_integral(e));
            out.points.push_back(po);
        }
    }

    // Planes
    {
        auto vplanes = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                                cardillo::PhysicsSystem::C_Position3,
                                cardillo::PhysicsSystem::C_Plane>();
        for (auto [e, pos, pl] : vplanes.each()) {
            PlaneOut po; po.center = pos.value; po.normal = pl.normal; po.up = pl.up; po.sizeX = pl.sizeX; po.sizeY = pl.sizeY;
            out.planes.push_back(po);
        }
    }

    // Cubes
    {
        auto vcubes = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                               cardillo::PhysicsSystem::C_Position3,
                               cardillo::PhysicsSystem::C_Cube,
                               cardillo::PhysicsSystem::C_Orientation>();
        for (auto [e, pos, cb, ori] : vcubes.each()) {
            CubeOut co;
            Quaternion4r q_local = cb.q;
            Quaternion4r q_world = ori.value * q_local;
            co.center = pos.value + q_world * cb.center; // apply local center offset under local cube rotation
            co.halfExtents = cb.halfExtents; co.q = q_world;
            co.entityId = static_cast<int>(entt::to_integral(e));
            co.partition = partitionFromBodyIndex_(sys, reg, e);
            co.isDynamic = reg.any_of<cardillo::PhysicsSystem::C_PhysicsObject>(e);
            if (reg.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(e) && reg.any_of<cardillo::PhysicsSystem::C_AngularVelocity3>(e)) {
                const auto& vlin = reg.get<cardillo::PhysicsSystem::C_LinearVelocity3>(e).value;
                const auto& omega_body = reg.get<cardillo::PhysicsSystem::C_AngularVelocity3>(e).value;
                co.vlin = vlin; co.omega = omega_body; co.hasKinematics = true;
            }
            out.cubes.push_back(co);
        }
    }

    // Capsules (emitted as triangle meshes)
    {
        auto vcaps = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                               cardillo::PhysicsSystem::C_Capsule,
                               cardillo::PhysicsSystem::C_Position3,
                               cardillo::PhysicsSystem::C_Orientation>();
        for (auto [e, cap, pos, ori] : vcaps.each()) {
            std::vector<Vector3r> verts;
            std::vector<Eigen::Vector3i> tris;
            generateCapsuleMesh(8, 3, 1, cap.radius, cap.halfLength, verts, tris);
            if (verts.empty() || tris.empty()) continue;
            MeshOut mo;
            mo.entityId = static_cast<int>(entt::to_integral(e));
            mo.partition = partitionFromBodyIndex_(sys, reg, e);
            mo.center = pos.value;
            mo.isDynamic = reg.any_of<cardillo::PhysicsSystem::C_PhysicsObject>(e);
            Quaternion4r qn = ori.value; qn.normalize();
            const Matrix33r R = qn.toRotationMatrix();
            mo.vertices.reserve(verts.size());
            for (const auto& v : verts) {
                mo.vertices.push_back(R * v + pos.value);
            }
            mo.triangles = std::move(tris);
            if (reg.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(e) && reg.any_of<cardillo::PhysicsSystem::C_AngularVelocity3>(e)) {
                mo.vlin = reg.get<cardillo::PhysicsSystem::C_LinearVelocity3>(e).value;
                const auto& omega_body = reg.get<cardillo::PhysicsSystem::C_AngularVelocity3>(e).value;
                mo.omega = omega_body;
                mo.hasKinematics = true;
            }
            out.meshes.push_back(std::move(mo));
        }
    }

    // Meshes
    {
        auto vmeshes = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                                 cardillo::PhysicsSystem::C_Mesh,
                                 cardillo::PhysicsSystem::C_Position3,
                                 cardillo::PhysicsSystem::C_Orientation>();
        for (auto [e, cm, pos, ori] : vmeshes.each()) {
            MeshOut mo;
            mo.entityId = static_cast<int>(entt::to_integral(e));
            mo.partition = partitionFromBodyIndex_(sys, reg, e);
            mo.center = pos.value;
            mo.isDynamic = reg.any_of<cardillo::PhysicsSystem::C_PhysicsObject>(e);
            try {
                const auto& asset = sys.getMeshAsset(e);
                coal::BVHModelPtr_t bvh = asset.bvh;
                if (bvh && bvh->vertices && bvh->tri_indices) {
                    const auto& V = *bvh->vertices;
                    const auto& F = *bvh->tri_indices;
                    Quaternion4r qn = ori.value; qn.normalize();
                    const Matrix33r R = qn.toRotationMatrix();
                    mo.vertices.reserve(V.size());
                    for (const auto& v : V) {
                        Vector3r p((real_t)v[0], (real_t)v[1], (real_t)v[2]);
                        Vector3r pw = R * p + pos.value;
                        mo.vertices.push_back(pw);
                    }
                    mo.triangles.reserve(F.size());
                    for (const auto& t : F) mo.triangles.emplace_back((int)t[0], (int)t[1], (int)t[2]);
                    if (asset.hasUV && asset.uvs.size() == V.size()) { mo.uvs = asset.uvs; mo.hasUV = true; }
                    // Kinematics if available
                    if (reg.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(e) && reg.any_of<cardillo::PhysicsSystem::C_AngularVelocity3>(e)) {
                        const auto& vlin = reg.get<cardillo::PhysicsSystem::C_LinearVelocity3>(e).value;
                        const auto& omega_body = reg.get<cardillo::PhysicsSystem::C_AngularVelocity3>(e).value;
                        mo.vlin = vlin; mo.omega = omega_body; mo.hasKinematics = true; mo.R = R;
                    }
                    out.meshes.push_back(std::move(mo));
                }
            } catch (...) { /* skip */ }
        }
    }

    // Softbody surfaces (deformed meshes driven by point-mass nodes)
    {
        auto vsb = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                             cardillo::PhysicsSystem::C_SoftBodySurface>();
        for (auto [e, surf] : vsb.each()) {
            // Build a MeshOut by sampling current positions of the nodes
            MeshOut mo;
            mo.entityId = static_cast<int>(entt::to_integral(e));
            mo.partition = -1.0f; // visualization-only
            mo.center = Vector3r::Zero();
            mo.isDynamic = true; // nodes move over time
            // Gather vertex positions from node entities
            mo.vertices.reserve(surf.nodes.size());
            mo.perVertexVelocity.reserve(surf.nodes.size());
            const auto& reg2 = reg;
            for (entt::entity nodeEnt : surf.nodes) {
                if (reg2.valid(nodeEnt) && reg2.any_of<cardillo::PhysicsSystem::C_Position3>(nodeEnt)) {
                    mo.vertices.push_back(reg2.get<cardillo::PhysicsSystem::C_Position3>(nodeEnt).value);
                } else {
                    mo.vertices.emplace_back(Vector3r::Zero());
                }
                // Per-vertex velocity (default to zero if missing)
                if (reg2.valid(nodeEnt) && reg2.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(nodeEnt)) {
                    mo.perVertexVelocity.push_back(reg2.get<cardillo::PhysicsSystem::C_LinearVelocity3>(nodeEnt).value);
                } else {
                    mo.perVertexVelocity.emplace_back(Vector3r::Zero());
                }
            }
            mo.hasPerVertexVelocity = (mo.perVertexVelocity.size() == mo.vertices.size());
            // Copy triangles directly
            mo.triangles = surf.triangles;
            mo.hasUV = false;
            out.meshes.push_back(std::move(mo));
        }
    }

    // HeightFields (emit as meshes with UVs, decimated by stride)
    {
        auto vhfs = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                              cardillo::PhysicsSystem::C_HeightField,
                              cardillo::PhysicsSystem::C_HeightFieldVisualTag,
                              cardillo::PhysicsSystem::C_Position3,
                              cardillo::PhysicsSystem::C_Orientation>();
        for (auto [e, ch, pos, ori] : vhfs.each()) {
            try {
                const auto& asset = sys.getHeightFieldAsset(e);
                if (!asset.hf) continue;
                const auto& hf = *asset.hf;
                const auto& X = hf.getXGrid();
                const auto& Y = hf.getYGrid();
                const auto& H = hf.getHeights();
                const int rows = (int)H.rows();
                const int cols = (int)H.cols();
                if (rows <= 1 || cols <= 1) continue;
                const int s = std::max(1, m_hfStride);

                // Build index maps honoring stride and forcing inclusion of last indices
                std::vector<int> xIdx; xIdx.reserve((cols + s - 1) / s + 1);
                std::vector<int> yIdx; yIdx.reserve((rows + s - 1) / s + 1);
                for (int x = 0; x < cols; x += s) xIdx.push_back(x);
                if (xIdx.back() != cols - 1) xIdx.push_back(cols - 1);
                for (int y = 0; y < rows; y += s) yIdx.push_back(y);
                if (yIdx.back() != rows - 1) yIdx.push_back(rows - 1);
                const int nxs = (int)xIdx.size();
                const int nys = (int)yIdx.size();

                MeshOut mo;
                mo.entityId = static_cast<int>(entt::to_integral(e));
                mo.partition = -1.f;
                mo.center = pos.value;
                mo.isDynamic = false;

                Quaternion4r qn = ori.value; qn.normalize();
                const Matrix33r R = qn.toRotationMatrix();
                mo.vertices.reserve((size_t)nys * (size_t)nxs);
                mo.uvs.reserve((size_t)nys * (size_t)nxs);

                for (int yi = 0; yi < nys; ++yi) {
                    int y = yIdx[yi];
                    const float vf = nys > 1 ? (float)yi / (float)(nys - 1) : 0.f;
                    for (int xi = 0; xi < nxs; ++xi) {
                        int x = xIdx[xi];
                        const float uf = nxs > 1 ? (float)xi / (float)(nxs - 1) : 0.f;
                        Vector3r p_local((real_t)X[x], (real_t)Y[y], (real_t)H(y, x));
                        Vector3r p_world = R * p_local + pos.value;
                        mo.vertices.push_back(p_world);
                        mo.uvs.emplace_back(uf, 1.0 - vf);
                    }
                }

                // Triangulate grid: two triangles per strided cell
                mo.triangles.reserve((size_t)(nys - 1) * (size_t)(nxs - 1) * 2);
                for (int yi = 0; yi < nys - 1; ++yi) {
                    for (int xi = 0; xi < nxs - 1; ++xi) {
                        int i0 = yi * nxs + xi;
                        int i1 = i0 + 1;
                        int i2 = i0 + nxs;
                        int i3 = i2 + 1;
                        mo.triangles.emplace_back(i0, i1, i2);
                        mo.triangles.emplace_back(i1, i3, i2);
                    }
                }
                mo.hasUV = true;
                out.meshes.push_back(std::move(mo));
            } catch (...) {
                // ignore HF emit errors
            }
        }
    }

    // Rigid-body spheres as triangle meshes (UV-sphere)
    {
        // Pre-generate a shared unit sphere mesh
        static std::vector<Vector3r> unitVerts;
        static std::vector<Eigen::Vector3i> unitTris;
        static bool inited = false;
        if (!inited) { generateUVSphere(12, 24, unitVerts, unitTris); inited = true; }
        auto vspheres = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                                 cardillo::PhysicsSystem::C_RB_Sphere,
                                 cardillo::PhysicsSystem::C_Radius,
                                 cardillo::PhysicsSystem::C_Position3,
                                 cardillo::PhysicsSystem::C_Orientation>();
        for (auto [e, rad, pos, ori] : vspheres.each()) {
            MeshOut mo;
            mo.entityId = static_cast<int>(entt::to_integral(e));
            mo.partition = partitionFromBodyIndex_(sys, reg, e);
            mo.center = pos.value;
            mo.isDynamic = reg.any_of<cardillo::PhysicsSystem::C_PhysicsObject>(e);
            Quaternion4r qn = ori.value; qn.normalize();
            const Matrix33r R = qn.toRotationMatrix();
            const real_t r = rad.r;
            mo.vertices.reserve(unitVerts.size());
            for (const auto& v : unitVerts) {
                Vector3r p = R * (r * v) + pos.value;
                mo.vertices.push_back(p);
            }
            // Generate spherical UVs from the canonical unit sphere vertices
            mo.uvs.reserve(unitVerts.size());
            for (const auto& uvv : unitVerts) {
                const real_t vx = uvv.x();
                const real_t vy = uvv.y();
                const real_t vz = uvv.z();
                float u = static_cast<float>(std::atan2((double)vy, (double)vx) / (2.0 * M_PI) + 0.5);
                float v = static_cast<float>(std::acos((double)vz) / M_PI);
                mo.uvs.emplace_back(u, v);
            }
            mo.hasUV = true;
            mo.triangles = unitTris;
            // Kinematics
            if (reg.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(e) && reg.any_of<cardillo::PhysicsSystem::C_AngularVelocity3>(e)) {
                const auto& vlin = reg.get<cardillo::PhysicsSystem::C_LinearVelocity3>(e).value;
                const auto& omega_body = reg.get<cardillo::PhysicsSystem::C_AngularVelocity3>(e).value;
                mo.vlin = vlin; mo.omega = omega_body; mo.hasKinematics = true;
            }
            out.meshes.push_back(std::move(mo));
        }
    }

    return out;
}

// Big-endian helpers
inline uint32_t VtkWriterBinary::bswap32(uint32_t v) {
    return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8) | ((v & 0x00FF0000u) >> 8) | ((v & 0xFF000000u) >> 24);
}
inline void VtkWriterBinary::writeBE(std::ofstream& out, uint32_t v) {
    uint32_t b = bswap32(v); out.write(reinterpret_cast<const char*>(&b), sizeof(uint32_t));
}
inline void VtkWriterBinary::writeBE(std::ofstream& out, int32_t v) {
    writeBE(out, static_cast<uint32_t>(v));
}
inline void VtkWriterBinary::writeBE(std::ofstream& out, float v) {
    static_assert(sizeof(float) == 4, "float must be 4 bytes");
    uint32_t u; std::memcpy(&u, &v, 4); writeBE(out, u);
}

std::string VtkWriterBinary::buildPath(const std::string& prefix, int step) const {
    std::ostringstream ss; ss << prefix << '_' << std::setw(4) << std::setfill('0') << step << ".vtk";
    std::string filename = ss.str();
    return m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
}

void VtkWriterBinary::writeHeader(std::ofstream& out, const char* title) const {
    out << "# vtk DataFile Version 3.0\n";
    out << (title ? title : "Cardillo binary output") << "\n";
    out << "BINARY\n";
    out << "DATASET POLYDATA\n";
}

void VtkWriterBinary::writePointsBlock(std::ofstream& out, const Collected& data) const {
    std::size_t np = data.points.size();
    std::size_t nplanes = data.planes.size();
    std::size_t ncubes = data.cubes.size();
    std::size_t nmesh_pts = 0; for (const auto& m : data.meshes) nmesh_pts += m.vertices.size();
    const std::size_t ntotal = np + 4*nplanes + 8*ncubes + nmesh_pts;
    out << "POINTS " << ntotal << " float\n";

    // Points
    for (const auto& pt : data.points) { writeBE(out, f32(pt.pos.x())); writeBE(out, f32(pt.pos.y())); writeBE(out, f32(pt.pos.z())); }
    // Plane corners
    for (const auto& p : data.planes) {
        Vector3r n = p.normal; n.normalize();
        Vector3r a = p.up - p.up.dot(n)*n; if (a.norm() < 1e-12) a = Vector3r(1,0,0) - Vector3r(1,0,0).dot(n)*n; a.normalize();
        Vector3r b = n.cross(a);
        Vector3r c0 = p.center + (-p.sizeX)*a + (-p.sizeY)*b;
        Vector3r c1 = p.center + ( p.sizeX)*a + (-p.sizeY)*b;
        Vector3r c2 = p.center + ( p.sizeX)*a + ( p.sizeY)*b;
        Vector3r c3 = p.center + (-p.sizeX)*a + ( p.sizeY)*b;
        const Vector3r corners[4] = {c0,c1,c2,c3};
        for (int i = 0; i < 4; ++i) { writeBE(out, f32(corners[i].x())); writeBE(out, f32(corners[i].y())); writeBE(out, f32(corners[i].z())); }
    }
    // Cube corners
    for (const auto& c : data.cubes) {
        Matrix33r R = c.q.toRotationMatrix();
        Vector3r ex = R * Vector3r::UnitX();
        Vector3r ey = R * Vector3r::UnitY();
        Vector3r ez = R * Vector3r::UnitZ();
        Vector3r vx = ex * c.halfExtents.x();
        Vector3r vy = ey * c.halfExtents.y();
        Vector3r vz = ez * c.halfExtents.z();
        Vector3r p[8] = {
            c.center - vx - vy - vz,
            c.center + vx - vy - vz,
            c.center + vx + vy - vz,
            c.center - vx + vy - vz,
            c.center - vx - vy + vz,
            c.center + vx - vy + vz,
            c.center + vx + vy + vz,
            c.center - vx + vy + vz
        };
        for (int i = 0; i < 8; ++i) { writeBE(out, f32(p[i].x())); writeBE(out, f32(p[i].y())); writeBE(out, f32(p[i].z())); }
    }
    // Mesh vertices
    for (const auto& m : data.meshes) {
        for (const auto& p : m.vertices) { writeBE(out, f32(p.x())); writeBE(out, f32(p.y())); writeBE(out, f32(p.z())); }
    }
}

void VtkWriterBinary::writeVerticesBlock(std::ofstream& out, std::size_t np) const {
    out << "\nVERTICES " << np << ' ' << (2*np) << "\n";
    for (std::size_t i = 0; i < np; ++i) { writeBE(out, int32_t(1)); writeBE(out, int32_t(i)); }
}

void VtkWriterBinary::writePolygonsBlock(std::ofstream& out, const Collected& data) const {
    const std::size_t np = data.points.size();
    const std::size_t nplanes = data.planes.size();
    const std::size_t ncubes = data.cubes.size();
    std::size_t nmesh_tris = 0; for (const auto& m : data.meshes) nmesh_tris += m.triangles.size();

    const std::size_t nplane_pts = 4*nplanes;
    const std::size_t ncube_quads = 6*ncubes;
    const std::size_t nmesh_polys = nmesh_tris;
    const std::size_t npolys = nplanes + ncube_quads + nmesh_polys;
    const std::size_t polygonListSize = 5*(nplanes + ncube_quads) + 4*nmesh_polys;

    out << "POLYGONS " << npolys << ' ' << polygonListSize << "\n";

    // Planes as quads
    for (std::size_t i = 0; i < nplanes; ++i) {
        std::size_t base = np + 4*i;
        writeBE(out, int32_t(4));
        writeBE(out, int32_t(base+0)); writeBE(out, int32_t(base+1)); writeBE(out, int32_t(base+2)); writeBE(out, int32_t(base+3));
    }
    // Cubes as quads
    for (std::size_t i = 0; i < ncubes; ++i) {
        std::size_t base = np + nplane_pts + 8*i;
        const int faces[6][4] = {{0,1,2,3},{4,5,6,7},{0,3,7,4},{1,5,6,2},{0,4,5,1},{3,2,6,7}};
        for (int f = 0; f < 6; ++f) {
            writeBE(out, int32_t(4));
            for (int k = 0; k < 4; ++k) writeBE(out, int32_t(base + faces[f][k]));
        }
    }
    // Mesh triangles
    std::size_t base = np + nplane_pts + 8*ncubes;
    for (const auto& m : data.meshes) {
        for (const auto& tri : m.triangles) {
            writeBE(out, int32_t(3));
            writeBE(out, int32_t(base + tri[0])); writeBE(out, int32_t(base + tri[1])); writeBE(out, int32_t(base + tri[2]));
        }
        base += m.vertices.size();
    }
}

void VtkWriterBinary::writePointDataPts(std::ofstream& out, const Collected& data) const {
    const std::size_t np = data.points.size();
    out << "\nPOINT_DATA " << np << "\n";

    // mass
    out << "SCALARS mass float 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) writeBE(out, pt.mass);

    // velocity
    out << "VECTORS velocity float\n";
    for (const auto& pt : data.points) { writeBE(out, f32(pt.vel.x())); writeBE(out, f32(pt.vel.y())); writeBE(out, f32(pt.vel.z())); }

    // entity velocity: emit the entity's linear velocity (for particles this is the same as the point velocity)
    out << "VECTORS entity_velocity float\n";
    for (const auto& pt : data.points) { writeBE(out, f32(pt.vel.x())); writeBE(out, f32(pt.vel.y())); writeBE(out, f32(pt.vel.z())); }

    // radius
    out << "SCALARS radius float 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) writeBE(out, pt.radius);

    // partition
    out << "SCALARS partition float 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) writeBE(out, pt.partition);

    // entity_id
    out << "SCALARS entity_id int 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) writeBE(out, int32_t(pt.entityId));
}

void VtkWriterBinary::writePointDataGeo(std::ofstream& out, const Collected& data) const {
    const std::size_t nplanes = data.planes.size();
    const std::size_t ncubes = data.cubes.size();
    std::size_t nmesh_pts = 0; for (const auto& m : data.meshes) nmesh_pts += m.vertices.size();
    const std::size_t ntotal = 4*nplanes + 8*ncubes + nmesh_pts;

    out << "\nPOINT_DATA " << ntotal << "\n";

    // velocity: planes zero, cubes computed, meshes: use per-vertex velocities if present, else rigid approx/zero
    out << "VECTORS velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& cu : data.cubes) {
        Matrix33r R = cu.q.toRotationMatrix();
        Vector3r ex = R * Vector3r::UnitX();
        Vector3r ey = R * Vector3r::UnitY();
        Vector3r ez = R * Vector3r::UnitZ();
        Vector3r vx = ex * cu.halfExtents.x();
        Vector3r vy = ey * cu.halfExtents.y();
        Vector3r vz = ez * cu.halfExtents.z();
        Vector3r p[8] = {
            cu.center - vx - vy - vz,
            cu.center + vx - vy - vz,
            cu.center + vx + vy - vz,
            cu.center - vx + vy - vz,
            cu.center - vx - vy + vz,
            cu.center + vx - vy + vz,
            cu.center + vx + vy + vz,
            cu.center - vx + vy + vz
        };
        for (int i = 0; i < 8; ++i) {
            Vector3r v = Vector3r::Zero();
            if (cu.hasKinematics) { Vector3r r = p[i] - cu.center; v = cu.vlin + cu.omega.cross(r); }
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    }
    // Mesh vertices velocities: softbodies provide per-vertex values; otherwise use rigid kinematics if available
    for (const auto& m : data.meshes) {
        const bool usePV = m.hasPerVertexVelocity && (m.perVertexVelocity.size() == m.vertices.size());
        if (usePV) {
            for (const auto& v : m.perVertexVelocity) { writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z())); }
        } else {
            for (const auto& pw : m.vertices) {
                Vector3r v = Vector3r::Zero();
                if (m.hasKinematics) {
                    const Vector3r r_world = pw - m.center;
                    const Vector3r omega_world = m.R * m.omega;  // body → world
                    v = m.vlin + omega_world.cross(r_world);
                }
                writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
            }
        }
    }

    // angular velocity: planes zero, cubes repeated per-corner, meshes repeated per-vertex
    out << "VECTORS angular_velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& cu : data.cubes) {
        Vector3r omega = cu.hasKinematics ? cu.omega : Vector3r::Zero();
        for (int i = 0; i < 8; ++i) { writeBE(out, f32(omega.x())); writeBE(out, f32(omega.y())); writeBE(out, f32(omega.z())); }
    }
    for (const auto& m : data.meshes) {
        for (const auto& pw : m.vertices) {
            Vector3r omega = m.hasKinematics ? m.omega : Vector3r::Zero();
            writeBE(out, f32(omega.x())); writeBE(out, f32(omega.y())); writeBE(out, f32(omega.z()));
        }
    }

    // entity velocity: planes zero; cubes repeated per-corner; meshes: for softbodies reuse per-vertex velocity
    out << "VECTORS entity_velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& cu : data.cubes) {
        Vector3r vlin = cu.hasKinematics ? cu.vlin : Vector3r::Zero();
        for (int i = 0; i < 8; ++i) { writeBE(out, f32(vlin.x())); writeBE(out, f32(vlin.y())); writeBE(out, f32(vlin.z())); }
    }
    for (const auto& m : data.meshes) {
        const bool usePV = m.hasPerVertexVelocity && (m.perVertexVelocity.size() == m.vertices.size());
        if (usePV) {
            for (const auto& v : m.perVertexVelocity) { writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z())); }
        } else {
            for (const auto& pw : m.vertices) {
                Vector3r vlin = m.hasKinematics ? m.vlin : Vector3r::Zero();
                writeBE(out, f32(vlin.x())); writeBE(out, f32(vlin.y())); writeBE(out, f32(vlin.z()));
            }
        }
    }

    // partition: planes -1, cubes replicated, meshes replicated
    out << "SCALARS partition float 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) writeBE(out, -1.f);
    for (const auto& cu : data.cubes) { for (int i = 0; i < 8; ++i) writeBE(out, cu.partition); }
    for (const auto& m : data.meshes) { for (std::size_t i = 0; i < m.vertices.size(); ++i) writeBE(out, m.partition); }

    // entity_id
    out << "SCALARS entity_id int 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) writeBE(out, int32_t(-1));
    for (const auto& cu : data.cubes) { for (int i = 0; i < 8; ++i) writeBE(out, int32_t(cu.entityId)); }
    for (const auto& m : data.meshes) { for (std::size_t i = 0; i < m.vertices.size(); ++i) writeBE(out, int32_t(m.entityId)); }
}

// (Normals output removed intentionally)

void VtkWriterBinary::writeMeshTextureCoordinates(std::ofstream& out, const Collected& data) const {
    std::size_t nmesh_pts = 0; bool any = false;
    for (const auto& m : data.meshes) { nmesh_pts += m.vertices.size(); any = any || m.hasUV; }
    // Always emit a TCOORDS block for geometry so viewers can texture planes/cubes even if meshes lack UVs
    const std::size_t nplane_pts = 4*data.planes.size();
    const std::size_t ncube_pts = 8*data.cubes.size();
    out << "TEXTURE_COORDINATES tex 2 float\n";
    // Planes: emit UVs in meters, matching the corner order in writePointsBlock
    // Corner order is c0=(-hx,-hy), c1=(+hx,-hy), c2=(+hx,+hy), c3=(-hx,+hy) in the plane's local (a,b) frame.
    for (const auto& po : data.planes) {
        const float hx = static_cast<float>(po.sizeX);
        const float hy = static_cast<float>(po.sizeY);
        // Map local coordinates to [0, 2*hx] x [0, 2*hy]
        writeBE(out, 0.f);        writeBE(out, 0.f);        // c0: (-hx,-hy) -> (0,0)
        writeBE(out, 2.f*hx);     writeBE(out, 0.f);        // c1: (+hx,-hy) -> (2hx,0)
        writeBE(out, 2.f*hx);     writeBE(out, 2.f*hy);     // c2: (+hx,+hy) -> (2hx,2hy)
        writeBE(out, 0.f);        writeBE(out, 2.f*hy);     // c3: (-hx,+hy) -> (0,2hy)
    }

    // Cubes: emit simple corner UVs (u,v in {0,1}) matching the cube corner order used in writePointsBlock
    for (const auto& c : data.cubes) {
        // half extents in meters
        const float hx = static_cast<float>(c.halfExtents.x());
        const float hy = static_cast<float>(c.halfExtents.y());
        for (int i = 0; i < 8; ++i) {
            float lx = (i==1||i==2||i==5||i==6) ? hx : -hx;
            float ly = (i==2||i==3||i==6||i==7) ? hy : -hy;
            float u = lx + hx;
            float v = ly + hy;
            writeBE(out, u); writeBE(out, v);
        }
    }

    // Meshes: emit UVs or zeros
    for (const auto& m : data.meshes) {
        if (m.hasUV && m.uvs.size() == m.vertices.size()) {
            for (const auto& uv : m.uvs) { writeBE(out, uv.x()); writeBE(out, uv.y()); }
        } else {
            for (std::size_t i = 0; i < m.vertices.size(); ++i) { writeBE(out, 0.f); writeBE(out, 0.f); }
        }
    }

    // 
}

void VtkWriterBinary::writePointsOnly(int step, real_t /*time*/, const Collected& data) const {
    std::string path = buildPath(m_baseName + "_pts", step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;

    writeHeader(out, "Particles (binary)");
    // Only dynamic points: write their positions as full dataset (no planes/cubes/meshes)
    Collected onlyPts = data; onlyPts.planes.clear(); onlyPts.cubes.clear(); onlyPts.meshes.clear();
    writePointsBlock(out, onlyPts);
    writeVerticesBlock(out, onlyPts.points.size());
    writePointDataPts(out, onlyPts);
    out.close();
}

void VtkWriterBinary::writeGeometryOnly(int step, real_t /*time*/, const Collected& data) const {
    std::string path = buildPath(m_baseName + "_geo", step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;

    writeHeader(out, "Geometry (binary)");
    Collected onlyGeo; // no dynamic points
    onlyGeo.planes = data.planes;
    onlyGeo.cubes = data.cubes;
    onlyGeo.meshes = data.meshes;

    writePointsBlock(out, onlyGeo);
    writePolygonsBlock(out, onlyGeo);

    // POINT_DATA for planes/cubes/meshes
    writePointDataGeo(out, onlyGeo);
    writeMeshTextureCoordinates(out, onlyGeo);
    out.close();
}

void VtkWriterBinary::writeStaticGeometry(const Collected& data) const {
    // Write only static planes/cubes/meshes to a single file without step index
    std::string filename = m_baseName + std::string("_static_geo.vtk");
    std::string path = m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;
    writeHeader(out, "Static Geometry (binary)");
    Collected only;
    // Planes are static
    only.planes = data.planes;
    // Static cubes/meshes
    for (const auto& c : data.cubes) if (!c.isDynamic) only.cubes.push_back(c);
    for (const auto& m : data.meshes) if (!m.isDynamic) only.meshes.push_back(m);
    writePointsBlock(out, only);
    writePolygonsBlock(out, only);
    writePointDataGeo(out, only);
    writeMeshTextureCoordinates(out, only);
    out.close();
}

void VtkWriterBinary::writeDynamicGeometry(int step, const Collected& data) const {
    std::string path = buildPath(m_baseName + "_geo", step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;
    writeHeader(out, "Dynamic Geometry (binary)");
    Collected only;
    // No planes (static), include dynamic cubes/meshes only
    for (const auto& c : data.cubes) if (c.isDynamic) only.cubes.push_back(c);
    for (const auto& m : data.meshes) if (m.isDynamic) only.meshes.push_back(m);
    writePointsBlock(out, only);
    writePolygonsBlock(out, only);
    writePointDataGeo(out, only);
    writeMeshTextureCoordinates(out, only);
    out.close();
}

void VtkWriterBinary::writeStaticGeometryStep(int step, const Collected& data) const {
    // Mirror writeStaticGeometry but include step index so moving static objects are captured
    std::string path = buildPath(m_baseName + std::string("_static_geo"), step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;
    writeHeader(out, "Static Geometry (per-step) (binary)");
    Collected only;
    only.planes = data.planes; // planes treated as static
    for (const auto& c : data.cubes) if (!c.isDynamic) only.cubes.push_back(c);
    for (const auto& m : data.meshes) if (!m.isDynamic) only.meshes.push_back(m);
    writePointsBlock(out, only);
    writePolygonsBlock(out, only);
    writePointDataGeo(out, only);
    writeMeshTextureCoordinates(out, only);
    out.close();
}

void VtkWriterBinary::writeContacts(int step, const std::vector<cardillo::collision::Contact>& contacts, bool writeBodyVectors,
                                    cardillo::solver::WarmstartProvider* warmstartProvider) const {
    if (!m_writeContacts) return;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
    const std::string path = buildPath(m_contactsBase, step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;
    const std::size_t n = contacts.size();
    // Header and points
    writeHeader(out, "Collision contacts (binary)");
    out << "POINTS " << n << " float\n";
    for (const auto& c : contacts) { writeBE(out, f32(c.point.x())); writeBE(out, f32(c.point.y())); writeBE(out, f32(c.point.z())); }
    // Vertices section for points visibility
    out << "\nVERTICES " << n << ' ' << (2*n) << "\n";
    for (std::size_t i = 0; i < n; ++i) { writeBE(out, int32_t(1)); writeBE(out, int32_t(i)); }
    // Point data
    out << "\nPOINT_DATA " << n << "\n";
    out << "VECTORS normal float\n";
    for (const auto& c : contacts) { writeBE(out, f32(c.normal.x())); writeBE(out, f32(c.normal.y())); writeBE(out, f32(c.normal.z())); }
    out << "VECTORS tangent1 float\n";
    for (const auto& c : contacts) { writeBE(out, f32(c.tangent1.x())); writeBE(out, f32(c.tangent1.y())); writeBE(out, f32(c.tangent1.z())); }
    out << "VECTORS tangent2 float\n";
    for (const auto& c : contacts) { writeBE(out, f32(c.tangent2.x())); writeBE(out, f32(c.tangent2.y())); writeBE(out, f32(c.tangent2.z())); }
    if (writeBodyVectors) {
        out << "VECTORS normalA_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.normalA_body.x())); writeBE(out, f32(c.normalA_body.y())); writeBE(out, f32(c.normalA_body.z())); }
        out << "VECTORS normalB_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.normalB_body.x())); writeBE(out, f32(c.normalB_body.y())); writeBE(out, f32(c.normalB_body.z())); }
        out << "VECTORS pointA_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.pointA_body.x())); writeBE(out, f32(c.pointA_body.y())); writeBE(out, f32(c.pointA_body.z())); }
        out << "VECTORS pointB_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.pointB_body.x())); writeBE(out, f32(c.pointB_body.y())); writeBE(out, f32(c.pointB_body.z())); }
        out << "VECTORS tangent1A_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.tangent1A_body.x())); writeBE(out, f32(c.tangent1A_body.y())); writeBE(out, f32(c.tangent1A_body.z())); }
        out << "VECTORS tangent2A_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.tangent2A_body.x())); writeBE(out, f32(c.tangent2A_body.y())); writeBE(out, f32(c.tangent2A_body.z())); }
        out << "VECTORS tangent1B_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.tangent1B_body.x())); writeBE(out, f32(c.tangent1B_body.y())); writeBE(out, f32(c.tangent1B_body.z())); }
        out << "VECTORS tangent2B_body float\n";
        for (const auto& c : contacts) { writeBE(out, f32(c.tangent2B_body.x())); writeBE(out, f32(c.tangent2B_body.y())); writeBE(out, f32(c.tangent2B_body.z())); }
    }
    out << "SCALARS penetration float 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) writeBE(out, f32(c.penetration));

    // If a warmstart provider is available, emit percussion data (pn, tangential magnitude, and vector)
    if (warmstartProvider != nullptr) {
        // normal impulse magnitude
        out << "SCALARS pn float 1\nLOOKUP_TABLE default\n";
        for (const auto& c : contacts) {
            float pn = 0.f;
            if (c.prev_global_out_index >= 0) {
                auto opt = warmstartProvider->get(c.prev_global_out_index);
                if (opt) pn = static_cast<float>(opt->pn);
            }
            writeBE(out, pn);
        }

        // tangential magnitude (sqrt(pt1^2 + pt2^2))
        out << "SCALARS pt_mag float 1\nLOOKUP_TABLE default\n";
        for (const auto& c : contacts) {
            float ptmag = 0.f;
            if (c.prev_global_out_index >= 0) {
                auto opt = warmstartProvider->get(c.prev_global_out_index);
                if (opt) ptmag = static_cast<float>(std::sqrt((double)opt->pt1*opt->pt1 + (double)opt->pt2*opt->pt2));
            }
            writeBE(out, ptmag);
        }

        // percussion vector in world coordinates: pn*normal + pt1*tangent1 + pt2*tangent2
        out << "VECTORS percussion float\n";
        for (const auto& c : contacts) {
            Vector3r pvec = Vector3r::Zero();
            if (c.prev_global_out_index >= 0) {
                auto opt = warmstartProvider->get(c.prev_global_out_index);
                if (opt) {
                    pvec = (real_t)opt->pn * c.normal + (real_t)opt->pt1 * c.tangent1 + (real_t)opt->pt2 * c.tangent2;
                }
            }
            writeBE(out, f32(pvec.x())); writeBE(out, f32(pvec.y())); writeBE(out, f32(pvec.z()));
        }
    }
    out << "SCALARS id_a int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) writeBE(out, int32_t((int)entt::to_integral(c.a)));
    out << "SCALARS id_b int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) writeBE(out, int32_t((int)entt::to_integral(c.b)));
    out << "SCALARS friction_mu float 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) writeBE(out, f32(c.friction_mu));
    out << "SCALARS matched int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) writeBE(out, int32_t(c.prev_global_out_index >= 0 ? 1 : 0));
    out.close();
}

void VtkWriterBinary::writeSprings(int step, const cardillo::PhysicsSystem& sys) const {
    // Gather generic spring visuals (attachment B position, toAtoB) and
    // translation-rotation joint visuals (jointPos at B, frame, toAtoB)
    const auto& patterns = sys.constraintPatterns();
    std::vector<Vector3r> positions; positions.reserve(patterns.size());
    std::vector<Vector3r> toAtoB;   toAtoB.reserve(patterns.size());

    std::vector<Vector3r> tr_jointPos;
    std::vector<Vector3r> tr_toA;
    std::vector<Vector3r> tr_toB;
    std::vector<Vector3r> tr_ex;
    std::vector<Vector3r> tr_ey;
    std::vector<Vector3r> tr_ez;

    for (const auto& uptr : patterns) {
        if (!uptr) continue;

        // Generic spring: visualize at second attachment (xB),
        // and store vector from attachment B to attachment A.
        Vector3r xA, xB;
        if (uptr->getAttachPointsWorld(xA, xB)) {
            positions.push_back(xB);
            toAtoB.push_back(xA - xB);
        }

        // Translation-rotation style springs: extract joint frame from TranslationRotationConstraint
        if (auto* tr = dynamic_cast<const cardillo::physics::TranslationRotationConstraint*>(uptr.get())) {
            // Use the same attachment points: joint position is at xB,
            // and toA/toB are expressed from that point.
            const Vector3r jointPos = xB;

            // Build joint frame as concat of A's rotation and joint's A_K1J
            // A_IK1: world rotation of body A at current step
            const auto& reg = sys.ecs();
            const entt::entity a = tr->entityA();
            if (a == entt::null || !reg.all_of<cardillo::PhysicsSystem::C_Orientation>(a)) continue;
            const auto& qA = reg.get<cardillo::PhysicsSystem::C_Orientation>(a).value;
            const Matrix33r A_IK1 = qA.toRotationMatrix();

            const cardillo::physics::JointProperties& jp = tr->jointProperties();
            const Matrix33r A_IJ = A_IK1 * jp.A_K1J; // concat frame

            tr_jointPos.push_back(jointPos);
            tr_toA.push_back(xA - jointPos);
            tr_toB.push_back(xB - jointPos);
            tr_ex.push_back(A_IJ.col(0));
            tr_ey.push_back(A_IJ.col(1));
            tr_ez.push_back(A_IJ.col(2));
        }
    }

    const std::size_t n = positions.size();
    const std::size_t n_tr = tr_jointPos.size();
    if (n == 0 && n_tr == 0) return;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
    const std::string path = buildPath(m_springsBase, step);
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;
    // Header
    writeHeader(out, "Constraint springs (binary)");

    // Points: first generic springs (attachment B), then TR joint positions (also at B)
    const std::size_t totalPoints = n + n_tr;
    out << "POINTS " << totalPoints << " float\n";
    for (const auto& p : positions) { writeBE(out, f32(p.x())); writeBE(out, f32(p.y())); writeBE(out, f32(p.z())); }
    for (const auto& p : tr_jointPos) { writeBE(out, f32(p.x())); writeBE(out, f32(p.y())); writeBE(out, f32(p.z())); }

    // Vertices for visibility: one vertex per point
    out << "\nVERTICES " << totalPoints << ' ' << (2*totalPoints) << "\n";
    for (std::size_t i = 0; i < totalPoints; ++i) { writeBE(out, int32_t(1)); writeBE(out, int32_t(i)); }

    // Point data
    out << "\nPOINT_DATA " << totalPoints << "\n";

    // Generic toAtoB vectors (from attachment B to attachment A) for the
    // first n points; zeros for the translation-rotation-only points.
    out << "VECTORS toAtoB float\n";
    for (const auto& v : toAtoB) { writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z())); }
    for (std::size_t i = 0; i < n_tr; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }

    // Translation-rotation joint data for the last n_tr points; zeros for the first n
    auto writeVecField = [&](const char* name,
                             const std::vector<Vector3r>& prefixZeros,
                             const std::vector<Vector3r>& tail) {
        out << name << "\n";
        for (std::size_t i = 0; i < n; ++i) {
            writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f);
        }
        for (const auto& v : tail) {
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    };

    if (n_tr > 0) {
        writeVecField("VECTORS toA float", positions, tr_toA);
        writeVecField("VECTORS toB float", positions, tr_toB);
        writeVecField("VECTORS ex float",  positions, tr_ex);
        writeVecField("VECTORS ey float",  positions, tr_ey);
        writeVecField("VECTORS ez float",  positions, tr_ez);
    }
    out.close();
}

} // namespace cardillo::io
