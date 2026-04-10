#pragma once

#include <coal/BVH/BVH_model.h>
#include <coal/hfield.h>
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "../../misc/types.hpp"

namespace cardillo {

struct MeshAsset {
    coal::BVHModelPtr_t bvh;           // geometry BVH (possibly normalized)
    std::vector<Eigen::Vector2f> uvs;  // optional UVs parsed for OBJ
    bool hasUV = false;
    // Normalization metadata (unit density)
    Vector3r inertia_diag_unit = Vector3r::Zero();
    real_t volume = (real_t)0.0;            // signed volume
    Matrix33r Rpa = Matrix33r::Identity();  // principal-axes rotation
    Vector3r com = Vector3r::Zero();        // center of mass
    bool normalized = false;                // if true: COM-centered & PA-aligned
};

struct HeightFieldAsset {
    std::shared_ptr<coal::HeightField<coal::AABB>> hf;  // collider geometry
    int rows{0}, cols{0};
    real_t x_dim{1}, y_dim{1};
    real_t z_scale{1}, min_height{0};
    std::string path;
};

class PhysicsAssets {
   public:
    PhysicsAssets() = default;

    // Return mesh asset for path/scale. If normalized=true, asset is COM-centered & principal-axes
    // aligned.
    const MeshAsset& getMesh(const std::string& path, const Vector3r& scale = Vector3r::Ones(), bool normalized = false) const;

    const HeightFieldAsset& getHeightField(const std::string& exrPath, real_t x_dim, real_t y_dim, real_t z_scale = (real_t)1.0, real_t min_height = (real_t)0.0) const;

    void clear();

   private:
    mutable std::unordered_map<std::string, MeshAsset> m_meshCache;
    mutable std::unordered_map<std::string, HeightFieldAsset> m_hfCache;

    static std::string meshKey_(const std::string& path, const Vector3r& s, bool normalized);
    static std::string hfKey_(const std::string& path, real_t xd, real_t yd, real_t zs, real_t minh);

    // Builders
    MeshAsset buildMeshAsset_(const std::string& path, const Vector3r& scale, bool normalized) const;
    HeightFieldAsset buildHFAsset_(const std::string& path, real_t xd, real_t yd, real_t zs, real_t minh) const;
};

}  // namespace cardillo
