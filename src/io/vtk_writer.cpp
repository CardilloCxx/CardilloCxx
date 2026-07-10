#include "vtk_writer.hpp"
#include "../collision/collision_coal.hpp"
#include "../physics/assembly/dynamics_assembler.hpp"
#include "../physics/constraints/constraints.hpp"
#include "../rigid_body/rigid_body.hpp"
#include "mesh_generator.hpp"

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace cardillo::io {

namespace {

static inline float f32(real_t v) {
    return static_cast<float>(v);
}

static inline std::uint32_t entityKey(entt::entity e) {
    return static_cast<std::uint32_t>(entt::to_integral(e));
}

static real_t triangleArea(const Vector3r& a, const Vector3r& b, const Vector3r& c) {
    return static_cast<real_t>(0.5) * ((b - a).cross(c - a)).norm();
}

vtkSmartPointer<vtkFloatArray> makeFloatArray(const char* name, int numComponents) {
    vtkSmartPointer<vtkFloatArray> arr = vtkSmartPointer<vtkFloatArray>::New();
    arr->SetName(name);
    arr->SetNumberOfComponents(numComponents);
    return arr;
}

vtkSmartPointer<vtkIntArray> makeIntArray(const char* name) {
    vtkSmartPointer<vtkIntArray> arr = vtkSmartPointer<vtkIntArray>::New();
    arr->SetName(name);
    arr->SetNumberOfComponents(1);
    return arr;
}

}  // namespace

VtkWriter::VtkWriter(const cardillo::config::Config& cfg)
    : m_cfg(cfg),
      m_outputDir(cfg.output_folder),
      m_baseName(cfg.output_filename_prefix),
      m_frequency(cfg.output_interval_steps),
      m_writeSprings(true),
      m_springsBase(cfg.output_filename_prefix + std::string("_springs")),
      m_writeContactManifolds(cfg.output_write_contact_manifolds),
      m_contactManifoldsBase(cfg.output_filename_prefix + std::string("_contact_manifolds")) {
    if (m_frequency < 1) m_frequency = 1;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
    m_worker = std::thread(&VtkWriter::workerLoop, this);
}

VtkWriter::~VtkWriter() {
    {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        m_stop = true;
    }
    m_cvNotEmpty.notify_all();
    m_cvNotFull.notify_all();
    if (m_worker.joinable()) m_worker.join();
}

void VtkWriter::setOutputDir(const std::string& dir) {
    m_outputDir = dir;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
}

void VtkWriter::setBaseName(const std::string& name) {
    m_baseName = name;
}
void VtkWriter::setFrequency(int freq) {
    m_frequency = std::max(1, freq);
}

void VtkWriter::maybeWrite(int step, real_t time, const cardillo::World& sys, cardillo::collision::CollisionCoal* collision_mgr, cardillo::misc::TimingManager* timings,
                           cardillo::physics::DynamicsAssembler* dyn) {
    if (step % m_frequency == 0) write(step, time, sys, collision_mgr, timings, dyn);
}

void VtkWriter::write(int step, real_t time, const cardillo::World& sys, cardillo::collision::CollisionCoal* collision_mgr, cardillo::misc::TimingManager* timings,
                      cardillo::physics::DynamicsAssembler* dyn) {
    auto sc = timings->scope(cardillo::misc::TimingManager::TimerId::OutputWrite);

    std::vector<cardillo::collision::Contact> contacts;
    if (dyn) contacts = dyn->contacts();

    std::vector<EntityMesh> meshes;
    const auto& reg = sys.ecs();
    auto vis = reg.view<cardillo::C_VisualObject>();
    for (auto e : vis) {
        EntityMesh mesh;
        if (MeshGenerator::buildEntityMesh(sys, static_cast<entt::entity>(e), mesh)) {
            meshes.push_back(std::move(mesh));
        }
    }

    enrichPressure(meshes, sys, contacts);

    if (!m_staticGeoWritten || sys.isStructureDirty()) {
        std::vector<EntityMesh> onlyStatic;
        onlyStatic.reserve(meshes.size());
        for (const auto& m : meshes) {
            if (!m.isDynamic) onlyStatic.push_back(m);
        }
        const std::string filename = m_baseName + std::string("_static_geo.vtp");
        const std::string path = m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
        writeVtp(meshesToPolyData(onlyStatic), path);
        m_staticGeoWritten = true;
    }

    {
        std::vector<EntityMesh> onlyDynamic;
        onlyDynamic.reserve(meshes.size());
        for (const auto& m : meshes) {
            if (m.isDynamic) onlyDynamic.push_back(m);
        }
        const std::string path = buildPath(m_baseName + "_geo", step);
        enqueueFrame(meshesToPolyData(onlyDynamic), m_pvdGeo, pvdPath(m_baseName + "_geo"), path, fs::path(path).filename().string(), step, time);
    }

    if (m_writeContactManifolds && collision_mgr) {
        std::vector<cardillo::collision::ContactManifold> manifolds = collision_mgr->m_contactManifolds;
        enrichManifolds(manifolds, contacts, m_cfg.sim_dt);
        const std::string path = buildPath(m_contactManifoldsBase, step);
        enqueueFrame(contactManifoldsToPolyData(manifolds), m_pvdContactManifolds, pvdPath(m_contactManifoldsBase), path, fs::path(path).filename().string(), step, time);
    }

    if (m_writeSprings) {
        auto pd = springsToPolyData(sys);
        if (pd) {
            const std::string path = buildPath(m_springsBase, step);
            enqueueFrame(pd, m_pvdSprings, pvdPath(m_springsBase), path, fs::path(path).filename().string(), step, time);
        }
    }
}

void VtkWriter::enrichPressure(std::vector<EntityMesh>& meshes, const cardillo::World& sys, const std::vector<cardillo::collision::Contact>& contacts) const {
    const auto& reg = sys.ecs();
    const real_t invDt = (m_cfg.sim_dt > (real_t)0) ? ((real_t)1 / m_cfg.sim_dt) : (real_t)1;
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

void VtkWriter::enrichManifolds(std::vector<cardillo::collision::ContactManifold>& manifolds, const std::vector<cardillo::collision::Contact>& contacts, real_t dt) const {
    // detectAll() runs before the solve step, so Contact::last_impulse there is still last
    // step's warmstart seed. `contacts` here is fetched after the solve (same as the existing
    // contacts stream), so match manifold points to their solved Contact by position instead --
    // Contact::point is computed in makeContact() as the same (p1+p2)*0.5 mid-surface point that
    // coal::ContactPatch::getPoint() gives, so the match is exact (to floating-point noise).
    constexpr real_t kMatchDist2 = (real_t)1e-12;
    const real_t invDt = (dt > (real_t)0) ? ((real_t)1 / dt) : (real_t)1;

    for (auto& m : manifolds) {
        for (auto& pt : m.points) {
            for (const auto& c : contacts) {
                const bool samePair = (c.a == m.a && c.b == m.b) || (c.a == m.b && c.b == m.a);
                if (!samePair) continue;
                if ((c.point - pt.position).squaredNorm() <= kMatchDist2) {
                    pt.impulse = c.last_impulse;
                    break;
                }
            }
        }

        m.pressure = (real_t)0;
        if (m.points.size() < 3) continue;

        Vector3r areaVec = Vector3r::Zero();
        const std::size_t n = m.points.size();
        for (std::size_t i = 0; i < n; ++i) {
            areaVec += m.points[i].position.cross(m.points[(i + 1) % n].position);
        }
        const real_t area = (real_t)0.5 * std::abs(areaVec.dot(m.normal));
        if (area <= (real_t)1e-12) continue;

        real_t pnSum = (real_t)0;
        for (const auto& pt : m.points) pnSum += std::max<real_t>((real_t)0, pt.impulse.x());
        const real_t p = pnSum * invDt / area;
        m.pressure = std::isfinite(p) ? p : (real_t)0;
    }
}

std::string VtkWriter::buildPath(const std::string& prefix, int step) const {
    std::ostringstream ss;
    ss << prefix << '_' << std::setw(4) << std::setfill('0') << step << ".vtp";
    const std::string filename = ss.str();
    return m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
}

std::string VtkWriter::pvdPath(const std::string& prefix) const {
    const std::string filename = prefix + ".pvd";
    return m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
}

vtkSmartPointer<vtkPolyData> VtkWriter::meshesToPolyData(const std::vector<EntityMesh>& meshes) const {
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;

    auto velocity = makeFloatArray("velocity", 3);
    auto acceleration = makeFloatArray("acceleration", 3);
    auto angularVelocity = makeFloatArray("angular_velocity", 3);
    auto angularVelocityWorld = makeFloatArray("angular_velocity_world", 3);
    auto entityVelocity = makeFloatArray("entity_velocity", 3);
    auto entityId = makeIntArray("entity_id");
    auto beamLengthRatio = makeFloatArray("beam_length_ratio", 1);
    auto entityPressure = makeFloatArray("entity_pressure", 1);
    auto texCoords = makeFloatArray("tex", 2);

    std::size_t base = 0;
    for (const auto& m : meshes) {
        const bool usePV = m.perVertexVelocity.size() == m.vertices.size();
        const bool usePA = m.perVertexAcceleration.size() == m.vertices.size();
        const bool useUV = m.hasUV && m.uvs.size() == m.vertices.size();
        const Vector3r omegaWorld = m.R * m.omega;
        const Vector3r alphaWorld = m.R * m.alpha;

        for (std::size_t i = 0; i < m.vertices.size(); ++i) {
            const Vector3r& p = m.vertices[i];
            points->InsertNextPoint(f32(p.x()), f32(p.y()), f32(p.z()));

            Vector3r v = Vector3r::Zero();
            if (usePV) {
                v = m.perVertexVelocity[i];
            } else if (m.hasKinematics) {
                const Vector3r rWorld = p - m.center;
                v = m.vlin + omegaWorld.cross(rWorld);
            }
            velocity->InsertNextTuple3(f32(v.x()), f32(v.y()), f32(v.z()));

            Vector3r a = Vector3r::Zero();
            if (usePA) {
                a = m.perVertexAcceleration[i];
            } else if (m.hasKinematics) {
                const Vector3r rWorld = p - m.center;
                a = m.alin + alphaWorld.cross(rWorld) + omegaWorld.cross(omegaWorld.cross(rWorld));
            }
            acceleration->InsertNextTuple3(f32(a.x()), f32(a.y()), f32(a.z()));

            const Vector3r omega = m.hasKinematics ? m.omega : Vector3r::Zero();
            angularVelocity->InsertNextTuple3(f32(omega.x()), f32(omega.y()), f32(omega.z()));
            angularVelocityWorld->InsertNextTuple3(f32(omegaWorld.x()), f32(omegaWorld.y()), f32(omegaWorld.z()));

            const Vector3r ev = usePV ? m.perVertexVelocity[i] : (m.hasKinematics ? m.vlin : Vector3r::Zero());
            entityVelocity->InsertNextTuple3(f32(ev.x()), f32(ev.y()), f32(ev.z()));

            entityId->InsertNextValue(m.entityId);
            beamLengthRatio->InsertNextValue(f32(m.beamLengthRatio));
            entityPressure->InsertNextValue(f32(m.entityPressure));

            if (useUV) {
                texCoords->InsertNextTuple2(m.uvs[i].x(), m.uvs[i].y());
            } else {
                texCoords->InsertNextTuple2(0.0, 0.0);
            }
        }

        for (const auto& t : m.triangles) {
            if (t[0] < 0 || t[1] < 0 || t[2] < 0) continue;
            if ((std::size_t)t[0] >= m.vertices.size() || (std::size_t)t[1] >= m.vertices.size() || (std::size_t)t[2] >= m.vertices.size()) continue;
            vtkIdType tri[3] = {static_cast<vtkIdType>(base + (std::size_t)t[0]), static_cast<vtkIdType>(base + (std::size_t)t[1]), static_cast<vtkIdType>(base + (std::size_t)t[2])};
            polys->InsertNextCell(3, tri);
        }
        base += m.vertices.size();
    }

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(points);
    pd->SetPolys(polys);
    vtkPointData* pdData = pd->GetPointData();
    pdData->AddArray(velocity);
    pdData->AddArray(acceleration);
    pdData->AddArray(angularVelocity);
    pdData->AddArray(angularVelocityWorld);
    pdData->AddArray(entityVelocity);
    pdData->AddArray(entityId);
    pdData->AddArray(beamLengthRatio);
    pdData->AddArray(entityPressure);
    pdData->SetTCoords(texCoords);
    return pd;
}

vtkSmartPointer<vtkPolyData> VtkWriter::springsToPolyData(const cardillo::World& sys) const {
    const auto& patterns = sys.constraintPatterns();
    std::vector<Vector3r> positions;
    positions.reserve(patterns.size());
    std::vector<Vector3r> toAtoB;
    toAtoB.reserve(patterns.size());

    std::vector<Vector3r> tr_jointPos;
    std::vector<Vector3r> tr_toA;
    std::vector<Vector3r> tr_toB;
    std::vector<Vector3r> tr_ex;
    std::vector<Vector3r> tr_ey;
    std::vector<Vector3r> tr_ez;

    for (const auto& uptr : patterns) {
        if (!uptr) continue;

        Vector3r xA, xB;
        if (uptr->getAttachPointsWorld(xA, xB)) {
            positions.push_back(xB);
            toAtoB.push_back(xA - xB);
        }

        if (auto* tr = dynamic_cast<const cardillo::physics::TranslationRotationConstraint*>(uptr.get())) {
            const Vector3r jointPos = xB;

            const auto& reg = sys.ecs();
            const entt::entity a = tr->entityA();
            if (a == entt::null || !reg.all_of<cardillo::C_Orientation>(a)) continue;
            const auto stateA = cardillo::RigidBody::getState(reg, a);
            const Matrix33r A_IK1 = stateA.rotation;

            const cardillo::physics::JointProperties& jp = tr->jointProperties();
            const Matrix33r A_IJ = A_IK1 * jp.A_K1J;

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
    if (n == 0 && n_tr == 0) return nullptr;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> verts;
    const std::size_t totalPoints = n + n_tr;
    for (const auto& p : positions) points->InsertNextPoint(f32(p.x()), f32(p.y()), f32(p.z()));
    for (const auto& p : tr_jointPos) points->InsertNextPoint(f32(p.x()), f32(p.y()), f32(p.z()));
    for (std::size_t i = 0; i < totalPoints; ++i) {
        vtkIdType idx = static_cast<vtkIdType>(i);
        verts->InsertNextCell(1, &idx);
    }

    auto toAtoBArr = makeFloatArray("toAtoB", 3);
    for (const auto& v : toAtoB) toAtoBArr->InsertNextTuple3(f32(v.x()), f32(v.y()), f32(v.z()));
    for (std::size_t i = 0; i < n_tr; ++i) toAtoBArr->InsertNextTuple3(0.f, 0.f, 0.f);

    vtkSmartPointer<vtkFloatArray> toAArr, toBArr, exArr, eyArr, ezArr;
    if (n_tr > 0) {
        const auto fillTrField = [&](const char* name, const std::vector<Vector3r>& tail) {
            auto arr = makeFloatArray(name, 3);
            for (std::size_t i = 0; i < n; ++i) arr->InsertNextTuple3(0.f, 0.f, 0.f);
            for (const auto& v : tail) arr->InsertNextTuple3(f32(v.x()), f32(v.y()), f32(v.z()));
            return arr;
        };
        toAArr = fillTrField("toA", tr_toA);
        toBArr = fillTrField("toB", tr_toB);
        exArr = fillTrField("ex", tr_ex);
        eyArr = fillTrField("ey", tr_ey);
        ezArr = fillTrField("ez", tr_ez);
    }

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(points);
    pd->SetVerts(verts);
    vtkPointData* pdData = pd->GetPointData();
    pdData->AddArray(toAtoBArr);
    if (n_tr > 0) {
        pdData->AddArray(toAArr);
        pdData->AddArray(toBArr);
        pdData->AddArray(exArr);
        pdData->AddArray(eyArr);
        pdData->AddArray(ezArr);
    }
    return pd;
}

vtkSmartPointer<vtkPolyData> VtkWriter::contactManifoldsToPolyData(const std::vector<cardillo::collision::ContactManifold>& manifolds) const {
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> verts;
    vtkNew<vtkCellArray> lines;
    vtkNew<vtkCellArray> polys;

    auto normalArr = makeFloatArray("normal", 3);
    auto tangent1Arr = makeFloatArray("tangent1", 3);
    auto tangent2Arr = makeFloatArray("tangent2", 3);
    auto penetrationArr = makeFloatArray("penetration", 1);
    auto pressureArr = makeFloatArray("pressure", 1);
    auto frictionMuArr = makeFloatArray("friction_mu", 1);
    auto idAArr = makeIntArray("id_a");
    auto idBArr = makeIntArray("id_b");
    auto patchIdArr = makeIntArray("patch_id");
    auto pnArr = makeFloatArray("pn", 1);
    auto ptMagArr = makeFloatArray("pt_mag", 1);
    auto percussionArr = makeFloatArray("percussion", 3);

    vtkIdType base = 0;
    for (std::size_t im = 0; im < manifolds.size(); ++im) {
        const auto& m = manifolds[im];
        const std::size_t n = m.points.size();
        if (n == 0) continue;

        for (const auto& pt : m.points) {
            points->InsertNextPoint(f32(pt.position.x()), f32(pt.position.y()), f32(pt.position.z()));

            normalArr->InsertNextTuple3(f32(m.normal.x()), f32(m.normal.y()), f32(m.normal.z()));
            tangent1Arr->InsertNextTuple3(f32(m.tangent1.x()), f32(m.tangent1.y()), f32(m.tangent1.z()));
            tangent2Arr->InsertNextTuple3(f32(m.tangent2.x()), f32(m.tangent2.y()), f32(m.tangent2.z()));
            penetrationArr->InsertNextValue(f32(m.penetration));
            pressureArr->InsertNextValue(f32(m.pressure));
            frictionMuArr->InsertNextValue(f32(m.friction_mu));
            idAArr->InsertNextValue((int)entt::to_integral(m.a));
            idBArr->InsertNextValue((int)entt::to_integral(m.b));
            patchIdArr->InsertNextValue((int)im);

            const real_t pn = std::max<real_t>((real_t)0, pt.impulse.x());
            const real_t pt1 = pt.impulse.y();
            const real_t pt2 = pt.impulse.z();
            pnArr->InsertNextValue(f32(pn));
            ptMagArr->InsertNextValue(static_cast<float>(std::sqrt((double)pt1 * pt1 + (double)pt2 * pt2)));
            const Vector3r pvec = pn * m.normal + pt1 * m.tangent1 + pt2 * m.tangent2;
            percussionArr->InsertNextTuple3(f32(pvec.x()), f32(pvec.y()), f32(pvec.z()));
        }

        if (n == 1) {
            vtkIdType idx = base;
            verts->InsertNextCell(1, &idx);
        } else if (n == 2) {
            vtkIdType ids[2] = {base, base + 1};
            lines->InsertNextCell(2, ids);
        } else {
            std::vector<vtkIdType> ids(n);
            for (std::size_t i = 0; i < n; ++i) ids[i] = base + static_cast<vtkIdType>(i);
            polys->InsertNextCell(static_cast<vtkIdType>(n), ids.data());
        }
        base += static_cast<vtkIdType>(n);
    }

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(points);
    pd->SetVerts(verts);
    pd->SetLines(lines);
    pd->SetPolys(polys);
    vtkPointData* pdData = pd->GetPointData();
    pdData->AddArray(normalArr);
    pdData->AddArray(tangent1Arr);
    pdData->AddArray(tangent2Arr);
    pdData->AddArray(penetrationArr);
    pdData->AddArray(pressureArr);
    pdData->AddArray(frictionMuArr);
    pdData->AddArray(idAArr);
    pdData->AddArray(idBArr);
    pdData->AddArray(patchIdArr);
    pdData->AddArray(pnArr);
    pdData->AddArray(ptMagArr);
    pdData->AddArray(percussionArr);
    return pd;
}

void VtkWriter::writeVtp(const vtkSmartPointer<vtkPolyData>& pd, const std::string& path) {
    if (!pd) return;
    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(path.c_str());
    writer->SetInputData(pd);
    writer->SetDataModeToBinary();
    writer->SetCompressorTypeToNone();
    writer->Write();
}

void VtkWriter::enqueueFrame(vtkSmartPointer<vtkPolyData> pd, PvdWriter& pvd, const std::string& pvdFilePath, const std::string& vtpPath, const std::string& vtpBaseName, int step,
                            real_t time) {
    if (!pd) return;
    std::unique_lock<std::mutex> lock(m_queueMutex);
    m_cvNotFull.wait(lock, [this] { return m_queue.size() < kMaxQueueDepth || m_stop; });
    if (m_stop) return;
    m_queue.push_back([pd = std::move(pd), &pvd, pvdFilePath, vtpPath, vtpBaseName, step, time]() {
        writeVtp(pd, vtpPath);
        pvd.addEntry(step, time, vtpBaseName);
        pvd.write(pvdFilePath);
    });
    lock.unlock();
    m_cvNotEmpty.notify_one();
}

void VtkWriter::workerLoop() {
    while (true) {
        std::function<void()> job;
        {
            std::unique_lock<std::mutex> lock(m_queueMutex);
            m_cvNotEmpty.wait(lock, [this] { return !m_queue.empty() || m_stop; });
            if (m_queue.empty()) {
                if (m_stop) return;
                continue;
            }
            job = std::move(m_queue.front());
            m_queue.pop_front();
        }
        m_cvNotFull.notify_one();
        job();
    }
}

}  // namespace cardillo::io
