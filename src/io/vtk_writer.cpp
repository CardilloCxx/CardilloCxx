#include "vtk_writer.hpp"

#include <filesystem>
#include <fstream>
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

void VtkWriter::maybeWrite(int step, const cardillo::PhysicsSystem& sys)
{
    if (step % m_frequency == 0) write(step, sys);
}

void VtkWriter::write(int step, const cardillo::PhysicsSystem& sys)
{
    // format step as 4-digit with leading zeros
    std::ostringstream ss;
    ss << m_baseName << '_' << std::setw(4) << std::setfill('0') << step << ".vtk";
    std::string filename = ss.str();
    std::string path = m_outputDir.empty() ? filename : (fs::path(m_outputDir) / filename).string();

    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return;

    // Gather visuals from ECS and compute sizes
    Collected data = collect(sys);
    const std::size_t np = data.points.size();
    const std::size_t nplanes = data.planes.size();
    const std::size_t ncubes = data.cubes.size();
    const std::size_t nplane_pts = 4 * nplanes;
    const std::size_t ncube_pts = 8 * ncubes;
    const std::size_t ntotal = np + nplane_pts + ncube_pts;

    // Legacy VTK ASCII PolyData with particle points + plane corners
    writeHeader(out, ntotal);
    writePoints(out, data);

    // Vertices connectivity: one vertex per particle point only
    writeVertices(out, np);

    // Planes and cubes as quads (polygons)
    writePolygons(out, data);

    // POINT_DATA for all points: particles get mass/velocity, planes/cubes get zeros
    writePointData(out, data);
    out.close();
}

VtkWriter::Collected VtkWriter::collect(const cardillo::PhysicsSystem& sys) const {
    Collected out;
    const auto& reg = sys.ecs();
    // Points: visual points with mass and velocity
    auto vpoints = reg.view<cardillo::PhysicsSystem::C_VisualObject,
                            cardillo::PhysicsSystem::C_PointVisualTag,
                            cardillo::PhysicsSystem::C_Position3,
                            cardillo::PhysicsSystem::C_Mass,
                            cardillo::PhysicsSystem::C_LinearVelocity3>();
    for (auto [e, pos, m, v] : vpoints.each()) {
        (void)e;
        out.points.emplace_back(pos.value, std::make_pair(static_cast<float>(m.m), v.value));
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
        cardillo::PhysicsSystem::Cube c; c.center = pos.value; c.halfExtents = cb.halfExtents;
        out.cubes.push_back(c);
    }
    return out;
}

void VtkWriter::writeHeader(std::ofstream& out, std::size_t ntotal) const {
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
        const float cx = static_cast<float>(c.center.x());
        const float cy = static_cast<float>(c.center.y());
        const float cz = static_cast<float>(c.center.z());
        const float hx = static_cast<float>(c.halfExtents.x());
        const float hy = static_cast<float>(c.halfExtents.y());
        const float hz = static_cast<float>(c.halfExtents.z());
        const float pts[8][3] = {
            {cx - hx, cy - hy, cz - hz},
            {cx + hx, cy - hy, cz - hz},
            {cx + hx, cy + hy, cz - hz},
            {cx - hx, cy + hy, cz - hz},
            {cx - hx, cy - hy, cz + hz},
            {cx + hx, cy - hy, cz + hz},
            {cx + hx, cy + hy, cz + hz},
            {cx - hx, cy + hy, cz + hz}
        };
        for (int i = 0; i < 8; ++i) {
            out << pts[i][0] << ' ' << pts[i][1] << ' ' << pts[i][2] << '\n';
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
    for (const auto& pt : data.points) out << pt.second.first << '\n';
    for (std::size_t i = 0; i < nplane_pts; ++i) out << 0.0f << '\n';
    for (std::size_t i = 0; i < ncube_pts; ++i) out << 0.0f << '\n';
    out << "VECTORS velocity float\n";
    for (const auto& pt : data.points) {
        const auto& v = pt.second.second;
        out << static_cast<float>(v.x()) << ' '
            << static_cast<float>(v.y()) << ' '
            << static_cast<float>(v.z()) << '\n';
    }
    for (std::size_t i = 0; i < nplane_pts; ++i) out << "0 0 0\n";
    for (std::size_t i = 0; i < ncube_pts; ++i) out << "0 0 0\n";
}

} // namespace cardillo::io
