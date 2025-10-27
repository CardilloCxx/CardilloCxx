#include "vtk_writer_binary.hpp"
#include "../collision/collision_coal.hpp"
#include "vtk_sphere_util.hpp"
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
    Collected data = collect(sys);
    writePointsOnly(step, 0, data);
    // Write static geometry once or when structure changed
    if (!m_staticGeoWritten || sys.isStructureDirty()) {
        writeStaticGeometry(data);
        m_staticGeoWritten = true;
    }
    // Always write dynamic geometry per step
    writeDynamicGeometry(step, data);
    if (m_writeContacts) {
        try {
            auto& mgr = const_cast<cardillo::PhysicsSystem&>(sys).collisionManager();
            if (sys.consumeStructureDirty()) mgr.rebuild();
            mgr.applyTransforms();
            auto contacts = mgr.detectAll();
            const bool writeBody = sys.config().output_contacts_body_vectors;
            writeContacts(step, contacts, writeBody);
        } catch (...) {
            // skip if collision manager not available
        }
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
        int bodiesPerRank = (Nb + worldSize - 1) / worldSize;
        return (float)cardillo::partitioning::NaivePartitioner::ownerOf(b, worldSize, bodiesPerRank);
    }
    return -1.f;
}

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
            (void)e;
            PlaneOut po; po.p.center = pos.value; po.p.normal = pl.normal; po.p.up = pl.up; po.p.sizeX = pl.sizeX; po.p.sizeY = pl.sizeY;
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
            co.cube.center = pos.value; co.cube.halfExtents = cb.halfExtents; co.cube.q = ori.q;
            co.entityId = static_cast<int>(entt::to_integral(e));
            co.partition = partitionFromBodyIndex_(sys, reg, e);
            co.isDynamic = reg.any_of<cardillo::PhysicsSystem::C_PhysicsObject>(e);
            if (reg.any_of<cardillo::PhysicsSystem::C_LinearVelocity3>(e) && reg.any_of<cardillo::PhysicsSystem::C_AngularVelocity3>(e)) {
                const auto& vlin = reg.get<cardillo::PhysicsSystem::C_LinearVelocity3>(e).value;
                const auto& omega_body = reg.get<cardillo::PhysicsSystem::C_AngularVelocity3>(e).value;
                Matrix33r R = co.cube.q.toRotationMatrix();
                Vector3r omega_world = R * omega_body;
                co.vlin = vlin; co.omega = omega_world; co.hasKinematics = true;
            }
            out.cubes.push_back(co);
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
                    const Matrix33r R = ori.q.toRotationMatrix();
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
                        Vector3r omega_world = R * omega_body;
                        mo.vlin = vlin; mo.omega = omega_world; mo.hasKinematics = true;
                    }
                    out.meshes.push_back(std::move(mo));
                }
            } catch (...) { /* skip */ }
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

                const Matrix33r R = ori.q.toRotationMatrix();
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
            const Matrix33r R = ori.q.toRotationMatrix();
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
                Vector3r omega_world = R * omega_body;
                mo.vlin = vlin; mo.omega = omega_world; mo.hasKinematics = true;
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
    for (const auto& po : data.planes) {
        const auto& p = po.p;
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
    for (const auto& co : data.cubes) {
        const auto& c = co.cube;
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

    // velocity: planes zero, cubes computed, meshes zero (rigid approx optional)
    out << "VECTORS velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& co : data.cubes) {
        const auto& cu = co.cube;
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
            if (co.hasKinematics) { Vector3r r = p[i] - cu.center; v = co.vlin + co.omega.cross(r); }
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    }
    // Mesh vertices velocities based on rigid kinematics if available
    for (const auto& m : data.meshes) {
        for (const auto& pw : m.vertices) {
            Vector3r v = Vector3r::Zero();
            if (m.hasKinematics) { Vector3r r = pw - m.center; v = m.vlin + m.omega.cross(r); }
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    }

    // angular velocity: planes zero, cubes repeated per-corner, meshes repeated per-vertex
    out << "VECTORS angular_velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& co : data.cubes) {
        Vector3r omega = co.hasKinematics ? co.omega : Vector3r::Zero();
        for (int i = 0; i < 8; ++i) { writeBE(out, f32(omega.x())); writeBE(out, f32(omega.y())); writeBE(out, f32(omega.z())); }
    }
    for (const auto& m : data.meshes) {
        for (const auto& pw : m.vertices) {
            Vector3r omega = m.hasKinematics ? m.omega : Vector3r::Zero();
            writeBE(out, f32(omega.x())); writeBE(out, f32(omega.y())); writeBE(out, f32(omega.z()));
        }
    }

    // entity velocity: planes zero, cubes repeated per-corner, meshes repeated per-vertex
    out << "VECTORS entity_velocity float\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f); }
    for (const auto& co : data.cubes) {
        Vector3r vlin = co.hasKinematics ? co.vlin : Vector3r::Zero();
        for (int i = 0; i < 8; ++i) { writeBE(out, f32(vlin.x())); writeBE(out, f32(vlin.y())); writeBE(out, f32(vlin.z())); }
    }
    for (const auto& m : data.meshes) {
        for (const auto& pw : m.vertices) {
            Vector3r vlin = m.hasKinematics ? m.vlin : Vector3r::Zero();
            writeBE(out, f32(vlin.x())); writeBE(out, f32(vlin.y())); writeBE(out, f32(vlin.z()));
        }
    }

    // partition: planes -1, cubes replicated, meshes replicated
    out << "SCALARS partition float 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) writeBE(out, -1.f);
    for (const auto& co : data.cubes) { for (int i = 0; i < 8; ++i) writeBE(out, co.partition); }
    for (const auto& m : data.meshes) { for (std::size_t i = 0; i < m.vertices.size(); ++i) writeBE(out, m.partition); }

    // entity_id
    out << "SCALARS entity_id int 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < 4*nplanes; ++i) writeBE(out, int32_t(-1));
    for (const auto& co : data.cubes) { for (int i = 0; i < 8; ++i) writeBE(out, int32_t(co.entityId)); }
    for (const auto& m : data.meshes) { for (std::size_t i = 0; i < m.vertices.size(); ++i) writeBE(out, int32_t(m.entityId)); }
}

void VtkWriterBinary::writeMeshTextureCoordinates(std::ofstream& out, const Collected& data) const {
    std::size_t nmesh_pts = 0; bool any = false;
    for (const auto& m : data.meshes) { nmesh_pts += m.vertices.size(); any = any || m.hasUV; }
    // If there are cubes present, we will emit UVs for them (simple 0/1 corner UVs).
    if (!any && data.cubes.empty()) return;
    const std::size_t nplane_pts = 4*data.planes.size();
    const std::size_t ncube_pts = 8*data.cubes.size();
    out << "TEXTURE_COORDINATES tex 2 float\n";
    // Planes: no UVs (emit zeros)
    for (std::size_t i = 0; i < nplane_pts; ++i) { writeBE(out, 0.f); writeBE(out, 0.f); }

    // Cubes: emit simple corner UVs (u,v in {0,1}) matching the cube corner order used in writePointsBlock
    for (const auto& c : data.cubes) {
        // half extents in meters
        const float hx = static_cast<float>(c.cube.halfExtents.x());
        const float hy = static_cast<float>(c.cube.halfExtents.y());
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

void VtkWriterBinary::writeContacts(int step, const std::vector<cardillo::collision::Contact>& contacts, bool writeBodyVectors) const {
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

} // namespace cardillo::io
