#include "vtk_writer_binary.hpp"
#include "../collision/collision_coal.hpp"
#include "mesh_generator.hpp"
#include "../physics/solver/warmstart.hpp"
#include "../physics/constraints/constraints.hpp"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace cardillo::io {

namespace {

static inline float f32(real_t v) { return static_cast<float>(v); }

static inline std::uint32_t entityKey(entt::entity e) {
    return static_cast<std::uint32_t>(entt::to_integral(e));
}

static real_t triangleArea(const Vector3r& a, const Vector3r& b, const Vector3r& c) {
    return static_cast<real_t>(0.5) * ((b - a).cross(c - a)).norm();
}

} // namespace

void VtkWriterBinary::setOutputDir(const std::string& dir) {
    m_outputDir = dir;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
}

void VtkWriterBinary::setBaseName(const std::string& name) { m_baseName = name; }
void VtkWriterBinary::setFrequency(int freq) { m_frequency = std::max(1, freq); }

void VtkWriterBinary::maybeWrite(int step, real_t time, const cardillo::World& sys,
                                 cardillo::collision::CollisionCoal* collision_mgr,
                                 cardillo::misc::TimingManager* timings,
                                 cardillo::physics::DynamicsAssembler* dyn) {
    if (step % m_frequency == 0) write(step, time, sys, collision_mgr, timings, dyn);
}

void VtkWriterBinary::write(int step, real_t /*time*/, const cardillo::World& sys,
                            cardillo::collision::CollisionCoal* collision_mgr,
                            cardillo::misc::TimingManager* timings,
                            cardillo::physics::DynamicsAssembler* dyn) {
    auto sc = timings->scope(cardillo::misc::TimingManager::TimerId::OutputWrite);

    std::vector<cardillo::collision::Contact> contacts;
    if (dyn) contacts = dyn->contacts();

    std::vector<EntityMesh> meshes;
    const auto& reg = sys.ecs();
    auto vis = reg.view<cardillo::C_VisualObject>();
    for (auto e : vis) {
        EntityMesh mesh;
        if (MeshGenerator::buildEntityMesh(sys, static_cast<entt::entity>(e), m_hfStride, mesh)) {
            meshes.push_back(std::move(mesh));
        }
    }

    enrichPressure(meshes, sys, contacts);

    if (!m_staticGeoWritten || sys.isStructureDirty()) {
        writeStaticGeometry(meshes);
        m_staticGeoWritten = true;
    }

    writeDynamicGeometry(step, meshes);

    if (m_writeContacts && collision_mgr) {
        const bool writeBody = m_cfg.output_contacts_body_vectors;
        writeContacts(step, contacts, writeBody);
    }

    if (m_writeSprings) {
        writeSprings(step, sys);
    }
}

void VtkWriterBinary::enrichPressure(std::vector<EntityMesh>& meshes,
                                     const cardillo::World& sys,
                                     const std::vector<cardillo::collision::Contact>& contacts) const {
    const auto& reg = sys.ecs();
    const real_t invDt = (m_cfg.sim_dt > (real_t)0)
                       ? ((real_t)1 / m_cfg.sim_dt)
                       : (real_t)1;
    constexpr real_t kMinArea = (real_t)1e-12;

    std::unordered_map<std::uint32_t, real_t> compressiveForce;
    compressiveForce.reserve(contacts.size() * 2 + 1);

    for (const auto& c : contacts) {
        const real_t pn = std::max((real_t)0, c.last_impulse(0)) * invDt;
        if (pn <= (real_t)0) continue;
        compressiveForce[entityKey(c.a)] += pn;
        compressiveForce[entityKey(c.b)] += pn;
    }

    for (auto& m : meshes) {
        m.entityPressure = (real_t)0;
        if (!m.isDynamic) continue;
        if (m.vertices.empty()) continue;
        if (!reg.valid(m.entity)) continue;

        real_t area = (real_t)0;
        for (const auto& t : m.triangles) {
            const int i0 = t[0];
            const int i1 = t[1];
            const int i2 = t[2];
            if (i0 < 0 || i1 < 0 || i2 < 0) continue;
            if ((std::size_t)i0 >= m.vertices.size() || (std::size_t)i1 >= m.vertices.size() || (std::size_t)i2 >= m.vertices.size()) continue;
            area += triangleArea(m.vertices[(std::size_t)i0], m.vertices[(std::size_t)i1], m.vertices[(std::size_t)i2]);
        }

        if (area <= kMinArea) continue;

        const auto it = compressiveForce.find(entityKey(m.entity));
        if (it == compressiveForce.end()) continue;

        const real_t p = it->second / area;
        m.entityPressure = std::isfinite(p) ? p : (real_t)0;
    }
}

inline uint32_t VtkWriterBinary::bswap32(uint32_t v) {
    return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8)
         | ((v & 0x00FF0000u) >> 8)  | ((v & 0xFF000000u) >> 24);
}

inline void VtkWriterBinary::writeBE(std::ofstream& out, uint32_t v) {
    uint32_t b = bswap32(v);
    out.write(reinterpret_cast<const char*>(&b), sizeof(uint32_t));
}

inline void VtkWriterBinary::writeBE(std::ofstream& out, int32_t v) {
    writeBE(out, static_cast<uint32_t>(v));
}

inline void VtkWriterBinary::writeBE(std::ofstream& out, float v) {
    static_assert(sizeof(float) == 4, "float must be 4 bytes");
    uint32_t u;
    std::memcpy(&u, &v, 4);
    writeBE(out, u);
}

std::string VtkWriterBinary::buildPath(const std::string& prefix, int step) const {
    std::ostringstream ss;
    ss << prefix << '_' << std::setw(4) << std::setfill('0') << step << ".vtk";
    const std::string filename = ss.str();
    return m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
}

void VtkWriterBinary::writeHeader(std::ofstream& out, const char* title) const {
    out << "# vtk DataFile Version 3.0\n";
    out << (title ? title : "Cardillo binary output") << "\n";
    out << "BINARY\n";
    out << "DATASET POLYDATA\n";
}

void VtkWriterBinary::writePointsBlock(std::ofstream& out, const std::vector<EntityMesh>& meshes) const {
    std::size_t nPoints = 0;
    for (const auto& m : meshes) nPoints += m.vertices.size();

    out << "POINTS " << nPoints << " float\n";
    for (const auto& m : meshes) {
        for (const auto& p : m.vertices) {
            writeBE(out, f32(p.x()));
            writeBE(out, f32(p.y()));
            writeBE(out, f32(p.z()));
        }
    }
}

void VtkWriterBinary::writePolygonsBlock(std::ofstream& out, const std::vector<EntityMesh>& meshes) const {
    std::size_t nPolys = 0;
    std::size_t listSize = 0;
    for (const auto& m : meshes) {
        nPolys += m.triangles.size();
        listSize += m.triangles.size() * 4;
    }

    out << "POLYGONS " << nPolys << ' ' << listSize << "\n";

    std::size_t base = 0;
    for (const auto& m : meshes) {
        for (const auto& t : m.triangles) {
            writeBE(out, int32_t(3));
            writeBE(out, int32_t(base + (std::size_t)t[0]));
            writeBE(out, int32_t(base + (std::size_t)t[1]));
            writeBE(out, int32_t(base + (std::size_t)t[2]));
        }
        base += m.vertices.size();
    }
}

void VtkWriterBinary::writePointDataGeo(std::ofstream& out, const std::vector<EntityMesh>& meshes) const {
    std::size_t ntotal = 0;
    for (const auto& m : meshes) ntotal += m.vertices.size();
    out << "\nPOINT_DATA " << ntotal << "\n";

    out << "VECTORS velocity float\n";
    for (const auto& m : meshes) {
        const bool usePV = m.perVertexVelocity.size() == m.vertices.size();
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            Vector3r v = Vector3r::Zero();
            if (usePV) {
                v = m.perVertexVelocity[i];
            } else if (m.hasKinematics) {
                const Vector3r rWorld = m.vertices[i] - m.center;
                const Vector3r omegaWorld = m.R * m.omega;
                v = m.vlin + omegaWorld.cross(rWorld);
            }
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    }

    out << "VECTORS acceleration float\n";
    for (const auto& m : meshes) {
        const bool usePA = m.perVertexAcceleration.size() == m.vertices.size();
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            Vector3r a = Vector3r::Zero();
            if (usePA) {
                a = m.perVertexAcceleration[i];
            } else if (m.hasKinematics) {
                const Vector3r omegaWorld = m.R * m.omega;
                const Vector3r alphaWorld = m.R * m.alpha;
                const Vector3r rWorld = m.vertices[i] - m.center;
                a = m.alin + alphaWorld.cross(rWorld) + omegaWorld.cross(omegaWorld.cross(rWorld));
            }
            writeBE(out, f32(a.x())); writeBE(out, f32(a.y())); writeBE(out, f32(a.z()));
        }
    }

    out << "VECTORS angular_velocity float\n";
    for (const auto& m : meshes) {
        const Vector3r omega = m.hasKinematics ? m.omega : Vector3r::Zero();
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            writeBE(out, f32(omega.x())); writeBE(out, f32(omega.y())); writeBE(out, f32(omega.z()));
        }
    }

    out << "VECTORS entity_velocity float\n";
    for (const auto& m : meshes) {
        const bool usePV = m.perVertexVelocity.size() == m.vertices.size();
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            Vector3r v = usePV ? m.perVertexVelocity[i] : (m.hasKinematics ? m.vlin : Vector3r::Zero());
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    }

    out << "SCALARS entity_id int 1\nLOOKUP_TABLE default\n";
    for (const auto& m : meshes) {
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            writeBE(out, int32_t(m.entityId));
        }
    }

    out << "SCALARS beam_length_ratio float 1\nLOOKUP_TABLE default\n";
    for (const auto& m : meshes) {
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            writeBE(out, f32(m.beamLengthRatio));
        }
    }

    out << "SCALARS entity_pressure float 1\nLOOKUP_TABLE default\n";
    for (const auto& m : meshes) {
        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            writeBE(out, f32(m.entityPressure));
        }
    }
}

void VtkWriterBinary::writeMeshTextureCoordinates(std::ofstream& out, const std::vector<EntityMesh>& meshes) const {
    out << "TEXTURE_COORDINATES tex 2 float\n";
    for (const auto& m : meshes) {
        const bool useUV = m.hasUV && m.uvs.size() == m.vertices.size();
        if (useUV) {
            for (const auto& uv : m.uvs) {
                writeBE(out, uv.x());
                writeBE(out, uv.y());
            }
        } else {
            for (std::size_t i = 0; i < m.vertices.size(); ++i) {
                writeBE(out, 0.f);
                writeBE(out, 0.f);
            }
        }
    }
}

void VtkWriterBinary::writeGeometryMeshList(const std::string& path,
                                            const std::vector<EntityMesh>& meshes,
                                            const char* title) const {
    std::ofstream out(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out) return;

    writeHeader(out, title);
    writePointsBlock(out, meshes);
    writePolygonsBlock(out, meshes);
    writePointDataGeo(out, meshes);
    writeMeshTextureCoordinates(out, meshes);
    out.close();
}

void VtkWriterBinary::writeStaticGeometry(const std::vector<EntityMesh>& meshes) const {
    std::vector<EntityMesh> only;
    only.reserve(meshes.size());
    for (const auto& m : meshes) {
        if (!m.isDynamic) only.push_back(m);
    }

    const std::string filename = m_baseName + std::string("_static_geo.vtk");
    const std::string path = m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
    writeGeometryMeshList(path, only, "Static Geometry (binary)");
}

void VtkWriterBinary::writeDynamicGeometry(int step, const std::vector<EntityMesh>& meshes) const {
    std::vector<EntityMesh> only;
    only.reserve(meshes.size());
    for (const auto& m : meshes) {
        if (m.isDynamic) only.push_back(m);
    }

    const std::string path = buildPath(m_baseName + "_geo", step);
    writeGeometryMeshList(path, only, "Dynamic Geometry (binary)");
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

    // Emit percussion data (pn, tangential magnitude, and vector) if any contact has a non-zero last_impulse
    bool hasImpulse = false;
    for (const auto &c : contacts) {
        if (c.last_impulse.squaredNorm() > (real_t)1e-12) { hasImpulse = true; break; }
    }
    if (hasImpulse) {
        // normal impulse magnitude
        out << "SCALARS pn float 1\nLOOKUP_TABLE default\n";
        for (const auto& c : contacts) {
            float pn = static_cast<float>(std::max<real_t>(c.last_impulse(0), (real_t)0));
            writeBE(out, pn);
        }

        // tangential magnitude (sqrt(pt1^2 + pt2^2))
        out << "SCALARS pt_mag float 1\nLOOKUP_TABLE default\n";
        for (const auto& c : contacts) {
            const real_t t1 = c.last_impulse(1);
            const real_t t2 = c.last_impulse(2);
            float ptmag = static_cast<float>(std::sqrt((double)t1 * t1 + (double)t2 * t2));
            writeBE(out, ptmag);
        }

        // percussion vector in world coordinates: pn*normal + pt1*tangent1 + pt2*tangent2
        out << "VECTORS percussion float\n";
        for (const auto& c : contacts) {
            Vector3r pvec = (real_t)c.last_impulse(0) * c.normal + (real_t)c.last_impulse(1) * c.tangent1 + (real_t)c.last_impulse(2) * c.tangent2;
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

void VtkWriterBinary::writeSprings(int step, const cardillo::World& sys) const {
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
            if (a == entt::null || !reg.all_of<cardillo::C_Orientation>(a)) continue;
            const auto& qA = reg.get<cardillo::C_Orientation>(a).value;
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
    auto writeVecField = [&](const char* name, const std::vector<Vector3r>& tail) {
        out << name << "\n";
        for (std::size_t i = 0; i < n; ++i) {
            writeBE(out, 0.f); writeBE(out, 0.f); writeBE(out, 0.f);
        }
        for (const auto& v : tail) {
            writeBE(out, f32(v.x())); writeBE(out, f32(v.y())); writeBE(out, f32(v.z()));
        }
    };

    if (n_tr > 0) {
        writeVecField("VECTORS toA float", tr_toA);
        writeVecField("VECTORS toB float", tr_toB);
        writeVecField("VECTORS ex float",  tr_ex);
        writeVecField("VECTORS ey float",  tr_ey);
        writeVecField("VECTORS ez float",  tr_ez);
    }
    out.close();
}

} // namespace cardillo::io
