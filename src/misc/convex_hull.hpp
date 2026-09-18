#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "types.hpp"

namespace cardillo::misc {

struct ConvexHull3D {
    struct Face {
        int v[3];
        Vector3r normal;
        real_t d;
        bool active{true};
    };

    static std::vector<Eigen::Vector3i> compute(const std::vector<Vector3r>& pts) {
        const int n = static_cast<int>(pts.size());
        if (n < 4) return {};

        // 1. Find 4 non-coplanar points to construct initial tetrahedron
        int p0 = 0, p1 = 0;
        real_t maxDistSq = 0;
        for (int i = 1; i < n; ++i) {
            real_t d = (pts[i] - pts[p0]).squaredNorm();
            if (d > maxDistSq) { maxDistSq = d; p1 = i; }
        }
        if (maxDistSq < static_cast<real_t>(1e-12)) return {};

        Vector3r lineDir = (pts[p1] - pts[p0]).normalized();
        int p2 = 0;
        real_t maxLineDist = 0;
        for (int i = 0; i < n; ++i) {
            if (i == p0 || i == p1) continue;
            real_t dist = ((pts[i] - pts[p0]).cross(lineDir)).squaredNorm();
            if (dist > maxLineDist) { maxLineDist = dist; p2 = i; }
        }
        if (maxLineDist < static_cast<real_t>(1e-12)) return {};

        Vector3r planeNormal = (pts[p1] - pts[p0]).cross(pts[p2] - pts[p0]).normalized();
        int p3 = 0;
        real_t maxPlaneDist = 0;
        for (int i = 0; i < n; ++i) {
            if (i == p0 || i == p1 || i == p2) continue;
            real_t dist = std::abs((pts[i] - pts[p0]).dot(planeNormal));
            if (dist > maxPlaneDist) { maxPlaneDist = dist; p3 = i; }
        }
        if (maxPlaneDist < static_cast<real_t>(1e-6)) return {}; // All points coplanar

        // Orient (p0, p1, p2) so p3 lies on the negative side
        if ((pts[p3] - pts[p0]).dot(planeNormal) > 0) {
            std::swap(p1, p2);
        }

        std::vector<Face> faces;
        auto addFace = [&](int v0, int v1, int v2) {
            Vector3r norm = (pts[v1] - pts[v0]).cross(pts[v2] - pts[v0]);
            real_t len = norm.norm();
            if (len < static_cast<real_t>(1e-12)) return;
            norm /= len;
            real_t d = -norm.dot(pts[v0]);
            faces.push_back(Face{{v0, v1, v2}, norm, d, true});
        };

        // Build initial tetrahedron faces (CCW outward)
        addFace(p0, p1, p2);
        addFace(p0, p2, p3);
        addFace(p0, p3, p1);
        addFace(p1, p3, p2);

        // 2. Incrementally add remaining points
        for (int i = 0; i < n; ++i) {
            if (i == p0 || i == p1 || i == p2 || i == p3) continue;

            const Vector3r& p = pts[i];
            
            // Find faces visible from point i
            std::vector<int> visibleIndices;
            for (size_t f = 0; f < faces.size(); ++f) {
                if (!faces[f].active) continue;
                if (faces[f].normal.dot(p) + faces[f].d > static_cast<real_t>(1e-7)) {
                    visibleIndices.push_back(static_cast<int>(f));
                }
            }

            if (visibleIndices.empty()) continue; // Point is inside or on boundary

            // Find horizon edges: edges shared by 1 visible face and 1 non-visible face
            struct Edge { int u, v; };
            std::vector<Edge> horizon;

            for (int fIdx : visibleIndices) {
                const auto& f = faces[fIdx];
                for (int e = 0; e < 3; ++e) {
                    int u = f.v[e];
                    int v = f.v[(e + 1) % 3];

                    bool neighborVisible = false;
                    for (int otherIdx : visibleIndices) {
                        if (otherIdx == fIdx) continue;
                        const auto& of = faces[otherIdx];
                        for (int oe = 0; oe < 3; ++oe) {
                            if (of.v[oe] == v && of.v[(oe + 1) % 3] == u) {
                                neighborVisible = true;
                                break;
                            }
                        }
                        if (neighborVisible) break;
                    }

                    if (!neighborVisible) {
                        horizon.push_back({u, v});
                    }
                }
            }

            // Deactivate visible faces
            for (int fIdx : visibleIndices) {
                faces[fIdx].active = false;
            }

            // Add new triangular faces connecting the horizon to the new point
            for (const auto& edge : horizon) {
                addFace(edge.u, edge.v, i);
            }
        }

        // Collect all active hull faces
        std::vector<Eigen::Vector3i> result;
        for (const auto& f : faces) {
            if (f.active) {
                result.emplace_back(f.v[0], f.v[1], f.v[2]);
            }
        }
        return result;
    }
};

}  // namespace
