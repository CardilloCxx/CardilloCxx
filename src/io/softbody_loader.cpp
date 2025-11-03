#include "softbody_loader.hpp"

#include <unordered_set>
#include <coal/mesh_loader/loader.h>
#include <coal/BVH/BVH_model.h>

namespace cardillo::io {

namespace {
inline uint64_t edge_key(int a, int b) {
    int i = a < b ? a : b;
    int j = a < b ? b : a;
    return (static_cast<uint64_t>(static_cast<uint32_t>(i)) << 32) |
           static_cast<uint32_t>(j);
}
}

bool load_obj_softbody(const std::string& path,
                       SoftBodyMesh& out,
                       const cardillo::Vector3r& scale) {
    out.positions.clear();
    out.edges.clear();
    out.triangles.clear();

    coal::MeshLoader loader(coal::BV_AABB);
    coal::BVHModelPtr_t bvh = loader.load(path);
    if (!bvh) return false;

    // Copy vertices with scaling
    if (bvh->vertices) {
        const auto& V = *bvh->vertices;
        out.positions.reserve(V.size());
        for (const auto& v : V) {
            out.positions.emplace_back(
                (real_t)v[0] * scale.x(),
                (real_t)v[1] * scale.y(),
                (real_t)v[2] * scale.z());
        }
    }

    // Build unique undirected edges from triangle indices if present, and record triangles
    if (bvh->tri_indices) {
        const auto& F = *bvh->tri_indices;
        std::unordered_set<uint64_t> Eset;
        Eset.reserve(F.size() * 3);
        for (const auto& tri : F) {
            const int a = tri[0];
            const int b = tri[1];
            const int c = tri[2];
            if (a >= 0 && b >= 0) Eset.insert(edge_key(a,b));
            if (b >= 0 && c >= 0) Eset.insert(edge_key(b,c));
            if (c >= 0 && a >= 0) Eset.insert(edge_key(c,a));
            // triangles
            if (a >= 0 && b >= 0 && c >= 0) {
                out.triangles.emplace_back(a, b, c);
            }
        }
        out.edges.reserve(Eset.size());
        for (uint64_t k : Eset) {
            int i = static_cast<int>(k >> 32);
            int j = static_cast<int>(k & 0xffffffffu);
            out.edges.emplace_back(i, j);
        }
    }
    return !out.positions.empty();
}

} // namespace cardillo::io
