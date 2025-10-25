#include "vtk_sphere_util.hpp"
#include <cmath>

namespace cardillo::io {

void generateUVSphere(int latSegments,
                      int lonSegments,
                      std::vector<Vector3r>& outVertices,
                      std::vector<Eigen::Vector3i>& outTriangles)
{
    if (latSegments < 3) latSegments = 3;
    if (lonSegments < 3) lonSegments = 3;
    outVertices.clear();
    outTriangles.clear();
    // Top pole
    outVertices.emplace_back(0,0,1);
    // Rings
    for (int i = 1; i < latSegments; ++i) {
        const real_t phi = (real_t)M_PI * (real_t)i / (real_t)latSegments; // (0,pi)
        const real_t z = std::cos(phi);
        const real_t r = std::sin(phi);
        for (int j = 0; j < lonSegments; ++j) {
            const real_t theta = (real_t)2.0 * (real_t)M_PI * (real_t)j / (real_t)lonSegments;
            const real_t x = r * std::cos(theta);
            const real_t y = r * std::sin(theta);
            outVertices.emplace_back(x, y, z);
        }
    }
    // Bottom pole
    const int bottomIndex = (int)outVertices.size();
    outVertices.emplace_back(0,0,-1);

    // Helper to index ring vertices
    auto ringIndex = [&](int ring, int seg) -> int {
        // ring in [0, latSegments-2] for interior rings
        // seg in [0, lonSegments-1]
        return 1 + ring * lonSegments + (seg % lonSegments);
    };

    // Top cap triangles
    for (int j = 0; j < lonSegments; ++j) {
        int a = 0; // top pole
        int b = ringIndex(0, j);
        int c = ringIndex(0, j+1);
        outTriangles.emplace_back(a, b, c);
    }
    // Middle quads split into triangles
    for (int i = 0; i < latSegments-2; ++i) {
        for (int j = 0; j < lonSegments; ++j) {
            int v00 = ringIndex(i, j);
            int v01 = ringIndex(i, j+1);
            int v10 = ringIndex(i+1, j);
            int v11 = ringIndex(i+1, j+1);
            outTriangles.emplace_back(v00, v10, v11);
            outTriangles.emplace_back(v00, v11, v01);
        }
    }
    // Bottom cap triangles
    const int lastRing = latSegments - 2;
    for (int j = 0; j < lonSegments; ++j) {
        int a = ringIndex(lastRing, j);
        int b = bottomIndex;
        int c = ringIndex(lastRing, j+1);
        outTriangles.emplace_back(a, b, c);
    }
}

}
