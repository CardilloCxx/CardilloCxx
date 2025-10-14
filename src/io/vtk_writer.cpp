#include "vtk_writer.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace fs = std::filesystem;

namespace cardillo::io {

VtkWriter::VtkWriter(std::string outputDir, std::string baseName, int frequency)
    : m_outputDir(std::move(outputDir)), m_baseName(std::move(baseName)), m_frequency(frequency)
{
    if (m_frequency < 1) m_frequency = 1;
    // create output directory if missing
    if (!m_outputDir.empty()) {
        fs::create_directories(m_outputDir);
    }
}

void VtkWriter::setOutputDir(const std::string& dir)
{
    m_outputDir = dir;
    if (!m_outputDir.empty()) {
        fs::create_directories(m_outputDir);
    }
}

void VtkWriter::setBaseName(const std::string& name) { m_baseName = name; }

void VtkWriter::setFrequency(int freq) { m_frequency = std::max(1, freq); }

void VtkWriter::maybeWrite(int step, real_t time, const cardillo::PhysicsSystem& sys)
{
    if (step % m_frequency == 0) { write(step, time, sys); }
}

void VtkWriter::write(int step, real_t time, const cardillo::PhysicsSystem& sys)
{
    std::cout << "VtkWriter: writing step " << step << " at time " << time << std::endl;

    Collected data = collect(sys);
    writePointsOnly(step, time, data);
    writeGeometryOnly(step, time, data);

    if (m_writeContacts) {
        auto contacts = cardillo::collision::detectAll(sys);
        writeContacts(step, contacts);
    }
}

VtkWriter::Collected VtkWriter::collect(const cardillo::PhysicsSystem& sys) const {
    Collected out;
    const auto& reg = sys.ecs();
    // Points: visual points with mass and velocity
    auto vpoints = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                            cardillo::PhysicsSystem::C_PointVisualTag,
                            cardillo::PhysicsSystem::C_Position3,
                            cardillo::PhysicsSystem::C_Mass,
                            cardillo::PhysicsSystem::C_LinearVelocity3,
                            cardillo::PhysicsSystem::C_Radius>();
    for (auto [e, pos, m, v, r] : vpoints.each()) {
        (void)e;
        out.points.emplace_back(pos.value, std::make_tuple(static_cast<float>(m.m), v.value, static_cast<float>(r.r)));
    }
    // Planes: project ECS plane component to public Plane struct
    auto vplanes = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                            cardillo::PhysicsSystem::C_PlaneVisualTag,
                            cardillo::PhysicsSystem::C_Position3,
                            cardillo::PhysicsSystem::C_Plane>();
    for (auto [e, pos, pl] : vplanes.each()) {
        (void)e;
        cardillo::PhysicsSystem::Plane p; p.center = pos.value; p.normal = pl.normal; p.up = pl.up; p.sizeX = pl.sizeX; p.sizeY = pl.sizeY;
        out.planes.push_back(p);
    }
    // Cubes
    auto vcubes = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                           cardillo::PhysicsSystem::C_CubeVisualTag,
                           cardillo::PhysicsSystem::C_Position3,
                           cardillo::PhysicsSystem::C_Cube>();
    for (auto [e, pos, cb] : vcubes.each()) {
        (void)e;
        cardillo::PhysicsSystem::Cube c; c.center = pos.value; c.halfExtents = cb.halfExtents; c.R = cb.R;
        out.cubes.push_back(c);
    }
    return out;
}

void VtkWriter::writeHeader(std::ofstream& out, real_t time, std::size_t ntotal) const {
    out << "# vtk DataFile Version 3.0\n";
    out << "Particles + planes exported by cardillo\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << ntotal << " float\n";
    out.setf(std::ios::fixed); out << std::setprecision(6);
}

void VtkWriter::writePoints(std::ofstream& out, const Collected& data) const {
    // particle points first
    for (const auto& pt : data.points) {
        const auto& c = pt.first;
        out << static_cast<float>(c.x()) << ' '
            << static_cast<float>(c.y()) << ' '
            << static_cast<float>(c.z()) << '\n';
    }
    // plane points appended
    for (const auto& p : data.planes) {
        Vector3r n = p.normal; n.normalize();
        Vector3r a = p.up - p.up.dot(n) * n; if (a.norm() < 1e-12) a = Vector3r(1,0,0) - Vector3r(1,0,0).dot(n)*n; a.normalize();
        Vector3r b = n.cross(a);
        Vector3r c0 = p.center + (-p.sizeX)*a + (-p.sizeY)*b;
        Vector3r c1 = p.center + ( p.sizeX)*a + (-p.sizeY)*b;
        Vector3r c2 = p.center + ( p.sizeX)*a + ( p.sizeY)*b;
        Vector3r c3 = p.center + (-p.sizeX)*a + ( p.sizeY)*b;
        const Vector3r corners[4] = {c0,c1,c2,c3};
        for (int i = 0; i < 4; ++i) {
            out << static_cast<float>(corners[i].x()) << ' '
                << static_cast<float>(corners[i].y()) << ' '
                << static_cast<float>(corners[i].z()) << '\n';
        }
    }
    // cube points appended
    for (const auto& c : data.cubes) {
        // Rotate local axes
        Vector3r ex = c.R * Vector3r::UnitX();
        Vector3r ey = c.R * Vector3r::UnitY();
        Vector3r ez = c.R * Vector3r::UnitZ();
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
        for (int i = 0; i < 8; ++i) {
            out << static_cast<float>(p[i].x()) << ' '
                << static_cast<float>(p[i].y()) << ' '
                << static_cast<float>(p[i].z()) << '\n';
        }
    }
}

void VtkWriter::writeVertices(std::ofstream& out, std::size_t np) const {
    out << "VERTICES " << np << ' ' << 2 * np << "\n";
    for (std::size_t i = 0; i < np; ++i) {
        out << 1 << ' ' << i << '\n';
    }
}

void VtkWriter::writePolygons(std::ofstream& out, const Collected& data) const {
    const std::size_t np = data.points.size();
    const std::size_t nplanes = data.planes.size();
    const std::size_t ncubes = data.cubes.size();
    const std::size_t nplane_pts = 4 * nplanes;
    const std::size_t ncube_quads = 6 * ncubes;
    const std::size_t npolys = nplanes + ncube_quads;
    out << "POLYGONS " << npolys << ' ' << 5 * npolys << "\n";
    for (std::size_t i = 0; i < nplanes; ++i) {
        std::size_t base = np + 4 * i;
        out << 4 << ' ' << base << ' ' << (base+1) << ' ' << (base+2) << ' ' << (base+3) << '\n';
    }
    for (std::size_t i = 0; i < ncubes; ++i) {
        std::size_t base = np + nplane_pts + 8 * i;
        out << 4 << ' ' << base+0 << ' ' << base+1 << ' ' << base+2 << ' ' << base+3 << '\n';
        out << 4 << ' ' << base+4 << ' ' << base+5 << ' ' << base+6 << ' ' << base+7 << '\n';
        out << 4 << ' ' << base+0 << ' ' << base+3 << ' ' << base+7 << ' ' << base+4 << '\n';
        out << 4 << ' ' << base+1 << ' ' << base+5 << ' ' << base+6 << ' ' << base+2 << '\n';
        out << 4 << ' ' << base+0 << ' ' << base+4 << ' ' << base+5 << ' ' << base+1 << '\n';
        out << 4 << ' ' << base+3 << ' ' << base+2 << ' ' << base+6 << ' ' << base+7 << '\n';
    }
}

void VtkWriter::writePointData(std::ofstream& out, const Collected& data) const {
    const std::size_t np = data.points.size();
    const std::size_t nplanes = data.planes.size();
    const std::size_t ncubes = data.cubes.size();
    const std::size_t nplane_pts = 4 * nplanes;
    const std::size_t ncube_pts = 8 * ncubes;
    const std::size_t ntotal = np + nplane_pts + ncube_pts;
    out << "\nPOINT_DATA " << ntotal << "\n";
    out << "SCALARS mass float 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) out << std::get<0>(pt.second) << '\n';
    for (std::size_t i = 0; i < nplane_pts; ++i) out << 0.0f << '\n';
    for (std::size_t i = 0; i < ncube_pts; ++i) out << 0.0f << '\n';

    out << "VECTORS velocity float\n";
    for (const auto& pt : data.points) {
        const auto& v = std::get<1>(pt.second);
        out << static_cast<float>(v.x()) << ' '
            << static_cast<float>(v.y()) << ' '
            << static_cast<float>(v.z()) << '\n';
    }
    for (std::size_t i = 0; i < nplane_pts; ++i) out << "0 0 0\n";
    for (std::size_t i = 0; i < ncube_pts; ++i) out << "0 0 0\n";

    out << "SCALARS radius float 1\nLOOKUP_TABLE default\n";
    for (const auto& pt : data.points) out << std::get<2>(pt.second) << '\n';
    for (std::size_t i = 0; i < nplane_pts; ++i) out << 0.0f << '\n';
    for (std::size_t i = 0; i < ncube_pts; ++i) out << 0.0f << '\n';
}

auto VtkWriter::computeCounts(const Collected& data) const -> Counts {
    Counts c;
    c.np = data.points.size();
    c.nplanes = data.planes.size();
    c.ncubes = data.cubes.size();
    c.nplane_pts = 4 * c.nplanes;
    c.ncube_pts = 8 * c.ncubes;
    c.ntotal = c.np + c.nplane_pts + c.ncube_pts;
    return c;
}

std::string VtkWriter::buildPath(const std::string& prefix, int step) const {
    std::ostringstream ss;
    ss << prefix << '_' << std::setw(4) << std::setfill('0') << step << ".vtk";
    std::string filename = ss.str();
    return m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();
}

void VtkWriter::writePointsOnly(int step, real_t time, const Collected& data) const {
    // dynamic points only to allow glyphs without affecting static geometry
    std::string path = buildPath(m_baseName + "_pts", step);
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return;
    const std::size_t np = data.points.size();
    writeHeader(out, time, np);
    Collected onlyPts = data; onlyPts.planes.clear(); onlyPts.cubes.clear();
    writePoints(out, onlyPts);
    writeVertices(out, np);
    writePointData(out, onlyPts);
    out.close();
}

void VtkWriter::writeGeometryOnly(int step, real_t time, const Collected& data) const {
    // static geometry only (planes, cubes)
    std::string path = buildPath(m_baseName + "_geo", step);
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return;
    const auto c = computeCounts(data);
    const std::size_t ntotal = c.nplane_pts + c.ncube_pts;
    writeHeader(out, time, ntotal);
    Collected onlyGeo; // empty points
    onlyGeo.planes = data.planes;
    onlyGeo.cubes = data.cubes;
    writePoints(out, onlyGeo);
    writePolygons(out, onlyGeo);
    // no POINT_DATA section for static geometry
    out.close();
}

void VtkWriter::writeContacts(int step, const std::vector<cardillo::collision::Contact>& contacts) const {
    if (!m_writeContacts) return;
    if (!m_outputDir.empty()) fs::create_directories(m_outputDir);
    const std::string path = buildPath(m_contactsBase, step);
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return;
    const std::size_t n = contacts.size();
    out << "# vtk DataFile Version 3.0\n";
    out << "Collision contacts exported by cardillo\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << n << " float\n";
    out.setf(std::ios::fixed); out << std::setprecision(6);
    for (const auto& c : contacts) {
        out << static_cast<float>(c.point.x()) << ' '
            << static_cast<float>(c.point.y()) << ' '
            << static_cast<float>(c.point.z()) << '\n';
    }
    out << "\nPOINT_DATA " << n << "\n";
    out << "VECTORS normal float\n";
    for (const auto& c : contacts) {
        out << static_cast<float>(c.normal.x()) << ' '
            << static_cast<float>(c.normal.y()) << ' '
            << static_cast<float>(c.normal.z()) << '\n';
    }
    out << "SCALARS penetration float 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) out << static_cast<float>(c.penetration) << '\n';
    out << "SCALARS id_a int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) out << static_cast<int>(entt::to_integral(c.a)) << '\n';
    out << "SCALARS id_b int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : contacts) out << static_cast<int>(entt::to_integral(c.b)) << '\n';
    out.close();
}

} // namespace cardillo::io
