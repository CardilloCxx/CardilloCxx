#include "pvd_writer.hpp"

#include <fstream>

namespace cardillo::io {

void PvdWriter::addEntry(int step, real_t time, const std::string& fileName) {
    m_entries.push_back({step, time, fileName});
}

void PvdWriter::write(const std::string& path) const {
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) return;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (const auto& e : m_entries) {
        out << "    <DataSet timestep=\"" << e.time << "\" part=\"0\" file=\"" << e.fileName << "\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
}

}  // namespace cardillo::io
