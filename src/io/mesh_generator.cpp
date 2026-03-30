#include "mesh_generator.hpp"
#include <cmath>

namespace cardillo::io {

void MeshGenerator::generateUVSphere(int latSegments,
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


void MeshGenerator::generateCapsuleMesh(int segmentsCircumference,
                         int segmentsHemisphere,
                         int segmentsCylinder,
                         real_t radius,
                         real_t halfLength,
                         std::vector<Vector3r>& vertices,
                         std::vector<Eigen::Vector3i>& triangles)
{
    vertices.clear();
    triangles.clear();
    if (radius <= (real_t)0) return;
    const int nCirc = std::max(segmentsCircumference, 3);
    const int nHem = std::max(segmentsHemisphere, 1);
    const int nCyl = std::max(segmentsCylinder, 1);

    const real_t pi = std::acos(static_cast<real_t>(-1));
    const real_t halfPi = pi * (real_t)0.5;

    auto addRing = [&](real_t z, real_t ringRadius) {
        std::vector<int> idx;
        idx.reserve(nCirc);
        for (int i = 0; i < nCirc; ++i) {
            const real_t angle = static_cast<real_t>(2.0 * pi * i / nCirc);
            const real_t cx = ringRadius * std::cos(angle);
            const real_t cy = ringRadius * std::sin(angle);
            idx.push_back(static_cast<int>(vertices.size()));
            vertices.emplace_back(cx, cy, z);
        }
        return idx;
    };

    std::vector<std::vector<int>> rings;
    rings.reserve(nHem * 2 + nCyl + 2);

    const int bottomPole = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, -halfLength - radius);

    for (int i = 1; i < nHem; ++i) {
        const real_t theta = (static_cast<real_t>(i) / nHem) * halfPi;
        const real_t ringRadius = radius * std::sin(theta);
        const real_t z = -halfLength - radius * std::cos(theta);
        rings.push_back(addRing(z, ringRadius));
    }

    rings.push_back(addRing(-halfLength, radius));

    for (int i = 1; i < nCyl; ++i) {
        const real_t t = static_cast<real_t>(i) / nCyl;
        const real_t z = -halfLength + t * (halfLength * 2);
        rings.push_back(addRing(z, radius));
    }

    rings.push_back(addRing(halfLength, radius));

    for (int i = nHem - 1; i >= 1; --i) {
        const real_t theta = (static_cast<real_t>(i) / nHem) * halfPi;
        const real_t ringRadius = radius * std::sin(theta);
        const real_t z = halfLength + radius * std::cos(theta);
        rings.push_back(addRing(z, ringRadius));
    }

    const int topPole = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, halfLength + radius);

    if (rings.empty()) return;

    const auto& firstRing = rings.front();
    for (int j = 0; j < nCirc; ++j) {
        const int a = bottomPole;
        const int b = firstRing[j];
        const int c = firstRing[(j + 1) % nCirc];
        triangles.emplace_back(a, b, c);
    }

    for (std::size_t i = 0; i + 1 < rings.size(); ++i) {
        const auto& r0 = rings[i];
        const auto& r1 = rings[i + 1];
        for (int j = 0; j < nCirc; ++j) {
            const int jn = (j + 1) % nCirc;
            triangles.emplace_back(r0[j], r1[j], r0[jn]);
            triangles.emplace_back(r0[jn], r1[j], r1[jn]);
        }
    }

    const auto& lastRing = rings.back();
    for (int j = 0; j < nCirc; ++j) {
        const int a = lastRing[j];
        const int b = lastRing[(j + 1) % nCirc];
        const int c = topPole;
        triangles.emplace_back(a, b, c);
    }
}

void MeshGenerator::generateCylinderMesh(int segmentsCircumference,
                          real_t radius,
                          real_t halfLength,
                          std::vector<Vector3r>& vertices,
                          std::vector<Eigen::Vector3i>& triangles)
{
    vertices.clear();
    triangles.clear();
    if (radius <= (real_t)0) return;
    const int nCirc = std::max(segmentsCircumference, 3);

    const real_t pi = std::acos(static_cast<real_t>(-1));

    const int bottomCenterIdx = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, -halfLength);
    const int topCenterIdx = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, +halfLength);

    std::vector<int> bottomRing;
    std::vector<int> topRing;
    bottomRing.reserve(nCirc);
    topRing.reserve(nCirc);

    for (int i = 0; i < nCirc; ++i) {
        const real_t angle = static_cast<real_t>(2.0 * pi * i / nCirc);
        const real_t cx = radius * std::cos(angle);
        const real_t cy = radius * std::sin(angle);
        bottomRing.push_back(static_cast<int>(vertices.size()));
        vertices.emplace_back(cx, cy, -halfLength);
        topRing.push_back(static_cast<int>(vertices.size()));
        vertices.emplace_back(cx, cy, +halfLength);
    }

    // Side faces
    for (int i = 0; i < nCirc; ++i) {
        const int i0 = bottomRing[i];
        const int i1 = bottomRing[(i + 1) % nCirc];
        const int j0 = topRing[i];
        const int j1 = topRing[(i + 1) % nCirc];
        triangles.emplace_back(i0, i1, j1);
        triangles.emplace_back(i0, j1, j0);
    }

    // Bottom cap
    for (int i = 0; i < nCirc; ++i) {
        const int i0 = bottomRing[i];
        const int i1 = bottomRing[(i + 1) % nCirc];
        triangles.emplace_back(bottomCenterIdx, i1, i0);
    }

    // Top cap
    for (int i = 0; i < nCirc; ++i) {
        const int i0 = topRing[i];
        const int i1 = topRing[(i + 1) % nCirc];
        triangles.emplace_back(topCenterIdx, i0, i1);
    }
}

void MeshGenerator::generateConeMesh(int segmentsCircumference,
                      real_t radius,
                      real_t height,
                      std::vector<Vector3r>& vertices,
                      std::vector<Eigen::Vector3i>& triangles)
{
    vertices.clear();
    triangles.clear();
    if (radius <= (real_t)0 || height <= (real_t)0) return;
    const int nCirc = std::max(segmentsCircumference, 3);

    const real_t pi = std::acos(static_cast<real_t>(-1));
    const real_t halfH = height * (real_t)0.5;

    const int baseCenterIdx = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, -halfH);

    std::vector<int> ringIdx;
    ringIdx.reserve(nCirc);
    for (int i = 0; i < nCirc; ++i) {
        const real_t angle = static_cast<real_t>(2.0 * pi * i / nCirc);
        const real_t cx = radius * std::cos(angle);
        const real_t cy = radius * std::sin(angle);
        ringIdx.push_back(static_cast<int>(vertices.size()));
        vertices.emplace_back(cx, cy, -halfH);
    }

    const int apexIdx = static_cast<int>(vertices.size());
    vertices.emplace_back(0, 0, halfH);

    // Base triangles
    for (int i = 0; i < nCirc; ++i) {
        const int i0 = ringIdx[i];
        const int i1 = ringIdx[(i + 1) % nCirc];
        triangles.emplace_back(baseCenterIdx, i1, i0);
    }

    // Side triangles
    for (int i = 0; i < nCirc; ++i) {
        const int i0 = ringIdx[i];
        const int i1 = ringIdx[(i + 1) % nCirc];
        triangles.emplace_back(i0, i1, apexIdx);
    }
}

} // namespace cardillo::io
