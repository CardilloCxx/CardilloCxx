#pragma once

#include <string>
#include <vector>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"
// Optional: contacts writing support
#include "../collision/types.hpp"

namespace cardillo::io {
    
// File format: .vtk (POLYDATA with VERTICES and POINT_DATA optional)
class VtkWriter {
public:
    // outputDir: folder to write files to
    // baseName: prefix of filenames (e.g., "points"); files will be baseName_0000.vtk, etc.
    // frequency: write every 'frequency' steps; if 1, write each step
    VtkWriter(std::string outputDir = ".", std::string baseName = "points", int frequency = 1);

    void setOutputDir(const std::string& dir);
    void setBaseName(const std::string& name);
    void setFrequency(int freq);
    void enableContactsOutput(bool enable, const std::string& base = "contacts") { m_writeContacts = enable; m_contactsBase = base; }

    void maybeWrite(int step, real_t time, const cardillo::PhysicsSystem& sys);
    void write(int step, real_t time, const cardillo::PhysicsSystem& sys);
private:
    // Helpers split from write for clarity
    struct Collected {
        // pos, (mass, vel, radius, partition)
        // partition is float to allow NaN when unknown
        std::vector<std::pair<Vector3r, std::tuple<float, Vector3r, float, float>>> points;
        std::vector<cardillo::PhysicsSystem::Plane> planes;
        std::vector<cardillo::PhysicsSystem::Cube> cubes;
    };

    struct Counts {
        std::size_t np = 0;
        std::size_t nplanes = 0;
        std::size_t ncubes = 0;
        std::size_t nplane_pts = 0;
        std::size_t ncube_pts = 0;
        std::size_t ntotal = 0;
    };

    Collected collect(const cardillo::PhysicsSystem& sys) const;
    Counts computeCounts(const Collected& data) const;
    std::string buildPath(const std::string& prefix, int step) const;
    void writeHeader(std::ofstream& out, real_t time, std::size_t ntotal) const;
    void writePoints(std::ofstream& out, const Collected& data) const;
    void writeVertices(std::ofstream& out, std::size_t np) const;
    void writePolygons(std::ofstream& out, const Collected& data) const;
    void writePointData(std::ofstream& out, const Collected& data) const;
    
    void writeContacts(int step, const std::vector<cardillo::collision::Contact>& contacts) const;


    // Split-output helpers
    void writePointsOnly(int step, real_t time, const Collected& data) const;
    void writeGeometryOnly(int step, real_t time, const Collected& data) const;

    std::string m_outputDir;
    std::string m_baseName;
    int m_frequency;
    bool m_writeContacts = false;
    std::string m_contactsBase = "contacts";
};

} // namespace cardillo::io
