#pragma once

#include <string>
#include <vector>
#include <utility>
#include "../misc/types.hpp"
#include <Eigen/Core>

namespace cardillo::io {

// Minimal soft-body mesh representation extracted from a triangle mesh:
// - positions: one point-mass per OBJ vertex
// - edges: unique, undirected edges derived from triangle topology
struct SoftBodyMesh {
    std::vector<cardillo::Vector3r> positions;
    std::vector<std::pair<int,int>> edges; // indices into positions
    std::vector<Eigen::Vector3i> triangles; // triangle indices into positions
};

// Load an OBJ mesh and extract a vertex-edge graph suitable for a mass-spring system.
// Returns true on success and fills out with vertex positions and unique edges.
// Optional uniform scale can be applied to positions.
bool load_obj_softbody(const std::string& path,
                       SoftBodyMesh& out,
                       const cardillo::Vector3r& scale = cardillo::Vector3r(1,1,1));

}
