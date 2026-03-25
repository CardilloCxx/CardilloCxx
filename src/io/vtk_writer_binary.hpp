#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include "misc/types.hpp"
#include "../physics/world.hpp"
#include "vtk_sphere_util.hpp"

namespace cardillo { namespace collision { struct Contact; } }

namespace cardillo::io {

// A clean binary VTK (legacy) writer mirroring VtkWriter's functionality
// - Writes two files per step: <base>_pts_####.vtk and <base>_geo_####.vtk
// - Uses POLYDATA dataset with POINTS and POLYGONS
// - PointData includes mass/velocity/radius/partition/entity_id (for points-only)
//   and velocity/partition/entity_id (+optional mesh UVs) for geometry
class VtkWriterBinary {
public:
    VtkWriterBinary(std::string outputDir = ".", std::string baseName = "scene_bin", int frequency = 1);

    void setOutputDir(const std::string& dir);
    void setBaseName(const std::string& name);
    void setFrequency(int freq);
    void setHeightFieldStride(int s) { m_hfStride = s > 0 ? s : 1; }

    void maybeWrite(int step, real_t time, const cardillo::World& sys);
    void write(int step, real_t time, const cardillo::World& sys);

    // Optional contacts output
    void enableContactsOutput(bool enable, const std::string& baseName) { m_writeContacts = enable; m_contactsBase = baseName; }
    // Optional springs output: writes one POINT per spring (attachment A) and a VECTORS field pointing to attachment B
    void enableSpringsOutput(bool enable, const std::string& baseName) { m_writeSprings = enable; m_springsBase = baseName; }

private:
    struct PointOut {
        Vector3r pos{Vector3r::Zero()};
        float mass{0.0f};
        Vector3r vel{Vector3r::Zero()};
        float radius{0.0f};
        float partition{-1.0f};
        int entityId{-1};
    };
    struct PlaneOut { Vector3r center{Vector3r::Zero()}; Vector3r normal{0,0,1}; Vector3r up{0,1,0}; real_t sizeX{5}, sizeY{5}; };
    struct CubeOut {
        Vector3r center{Vector3r::Zero()}; Vector3r halfExtents{0.5,0.5,0.5}; Quaternion4r q{Quaternion4r::Identity()};
        Vector3r vlin{Vector3r::Zero()};
        Vector3r omega{Vector3r::Zero()};
        float partition{-1.0f};
        int entityId{-1};
        float beamLengthRatio{1.0f};
        bool hasKinematics{false};
        bool isDynamic{false};
    };
    struct MeshOut {
        std::vector<Vector3r> vertices;            // world-space
        std::vector<Eigen::Vector3i> triangles;    // indices local to this mesh block
        std::vector<Eigen::Vector2f> uvs;          // optional
        bool hasUV{false};
        // Optional per-vertex velocities (used for softbody visual surfaces)
        std::vector<Vector3r> perVertexVelocity;   // same size as vertices when present
        bool hasPerVertexVelocity{false};
        // Kinematics for per-vertex velocities (rigid assumption)
        Vector3r center{Vector3r::Zero()};
        Vector3r vlin{Vector3r::Zero()};
        Vector3r omega{Vector3r::Zero()};
        Matrix33r R{Matrix33r::Identity()};
        bool hasKinematics{false};
        float partition{-1.0f};
        int entityId{-1};
        float beamLengthRatio{1.0f};
        bool isDynamic{false};
    };

    struct Collected {
        std::vector<PointOut> points;
        std::vector<PlaneOut> planes;
        std::vector<CubeOut> cubes;
        std::vector<MeshOut> meshes;
    };

    // Helpers
    Collected collect(const cardillo::World& sys) const;

    // Binary (big-endian) write helpers
    static inline uint32_t bswap32(uint32_t v);
    static inline void writeBE(std::ofstream& out, uint32_t v);
    static inline void writeBE(std::ofstream& out, int32_t v);
    static inline void writeBE(std::ofstream& out, float v);

    std::string buildPath(const std::string& prefix, int step) const;

    // Blocks
    void writeHeader(std::ofstream& out, const char* title) const;
    void writePointsBlock(std::ofstream& out, const Collected& data) const;
    void writeVerticesBlock(std::ofstream& out, std::size_t np) const;
    void writePolygonsBlock(std::ofstream& out, const Collected& data) const;
    void writePointDataPts(std::ofstream& out, const Collected& data) const;
    void writePointDataGeo(std::ofstream& out, const Collected& data) const;
    void writeMeshTextureCoordinates(std::ofstream& out, const Collected& data) const;

    // Contacts (binary legacy VTK)
    void writeContacts(int step, const std::vector<cardillo::collision::Contact>& contacts, bool writeBodyVectors,
                       cardillo::solver::WarmstartProvider* warmstartProvider) const;

    void writePointsOnly(int step, real_t time, const Collected& data) const;
    void writeGeometryOnly(int step, real_t time, const Collected& data) const;
    void writeStaticGeometry(const Collected& data) const;
    void writeDynamicGeometry(int step, const Collected& data) const;
    void writeSprings(int step, const cardillo::World& sys) const;
    void writeStaticGeometryStep(int step, const Collected& data) const; // write static geometry each timestep (entity transforms may change)

    std::string m_outputDir;
    std::string m_baseName;
    int m_frequency{1};
    int m_hfStride{8}; // decimation stride for heightfield VTK output
    bool m_writeContacts{false};
    std::string m_contactsBase{"contacts"};
    bool m_writeSprings{false};
    std::string m_springsBase{"springs"};
    mutable bool m_staticGeoWritten{false};
};

} // namespace cardillo::io
