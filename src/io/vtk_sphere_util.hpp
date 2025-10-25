#pragma once

#include <vector>
#include "../misc/types.hpp"

namespace cardillo::io {

// Generate a unit UV-sphere triangulation.
// latSegments >= 3 (rings including poles are latSegments+1 vertices in vertical)
// lonSegments >= 3
// Outputs vertices on the unit sphere (radius=1) in body-local frame and triangle indices.
void generateUVSphere(int latSegments,
                      int lonSegments,
                      std::vector<Vector3r>& outVertices,
                      std::vector<Eigen::Vector3i>& outTriangles);

}
