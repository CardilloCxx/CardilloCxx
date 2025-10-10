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

    const auto& masses = sys.masses();
    const std::size_t n = masses.size();

    // Legacy VTK ASCII PolyData with points and vertices
    out << "# vtk DataFile Version 3.0\n";
    out << "Point masses exported by cardillo\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << n << " float\n";
    out.setf(std::ios::fixed); out << std::setprecision(6);
    for (const auto& pm : masses) {
        out << static_cast<float>(pm.x.x()) << ' '
            << static_cast<float>(pm.x.y()) << ' '
            << static_cast<float>(pm.x.z()) << '\n';
    }

    // Vertices connectivity: one vertex per point
    out << "VERTICES " << n << ' ' << 2 * n << "\n";
    for (std::size_t i = 0; i < n; ++i) {
        out << 1 << ' ' << i << '\n';
    }

    // Optional: point data, e.g., mass and velocity
    out << "\nPOINT_DATA " << n << "\n";
    out << "SCALARS mass float 1\nLOOKUP_TABLE default\n";
    for (const auto& pm : masses) {
        out << static_cast<float>(pm.m) << '\n';
    }
    out << "VECTORS velocity float\n";
    for (const auto& pm : masses) {
        out << static_cast<float>(pm.v.x()) << ' '
            << static_cast<float>(pm.v.y()) << ' '
            << static_cast<float>(pm.v.z()) << '\n';
    }
    out.close();
}

} // namespace cardillo::io
