#pragma once

#include <vector>
#include "misc/types.hpp"

namespace cardillo::io {

class MeshGenerator {
public:
    // Generate a unit UV-sphere triangulation.
    static void generateUVSphere(int latSegments,
                                 int lonSegments,
                                 std::vector<Vector3r>& outVertices,
                                 std::vector<Eigen::Vector3i>& outTriangles);

    static void generateCapsuleMesh(int segmentsCircumference,
                                    int segmentsHemisphere,
                                    int segmentsCylinder,
                                    real_t radius,
                                    real_t halfLength,
                                    std::vector<Vector3r>& vertices,
                                    std::vector<Eigen::Vector3i>& triangles);

    static void generateCylinderMesh(int segmentsCircumference,
                                     real_t radius,
                                     real_t halfLength,
                                     std::vector<Vector3r>& vertices,
                                     std::vector<Eigen::Vector3i>& triangles);

    static void generateConeMesh(int segmentsCircumference,
                                 real_t radius,
                                 real_t height,
                                 std::vector<Vector3r>& vertices,
                                 std::vector<Eigen::Vector3i>& triangles);
};

}
