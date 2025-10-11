#pragma once

#include <string>
#include <vector>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"

namespace cardillo::io {

// Writes point positions to legacy VTK ASCII files (one file per call)
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

    // Call per step; will write if (step % frequency == 0)
    // Exports all point masses in the PhysicsSystem
    void maybeWrite(int step, const cardillo::PhysicsSystem& sys);

    // Force a write for a specific step index
    void write(int step, const cardillo::PhysicsSystem& sys);

private:
    // Helpers split from write for clarity
    struct Collected {
        std::vector<std::pair<Vector3r, std::pair<float, Vector3r>>> points; // pos, (mass, vel)
        std::vector<cardillo::PhysicsSystem::Plane> planes;
        std::vector<cardillo::PhysicsSystem::Cube> cubes;
    };

    Collected collect(const cardillo::PhysicsSystem& sys) const;
    void writeHeader(std::ofstream& out, std::size_t ntotal) const;
    void writePoints(std::ofstream& out, const Collected& data) const;
    void writeVertices(std::ofstream& out, std::size_t np) const;
    void writePolygons(std::ofstream& out, const Collected& data) const;
    void writePointData(std::ofstream& out, const Collected& data) const;

    std::string m_outputDir;
    std::string m_baseName;
    int m_frequency;
};

} // namespace cardillo::io
