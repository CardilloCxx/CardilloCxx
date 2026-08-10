#pragma once

#ifdef NORMAL
#undef NORMAL
#endif

#include <filament/Engine.h>
#include <filament/IndexBuffer.h>
#include <filament/VertexBuffer.h>
#include <filament/RenderableManager.h>
#include <vector>

#include <physics/assets/assets.hpp>
#include "geometry/Geometry.h"


namespace frontend {

    enum class PrimitiveType {
        Cube,
        Sphere,
        Capsule,
        Cylinder,
        Cone,
    };

    class Assets {
    public:
        Assets(filament::Engine& engine);
        ~Assets();
        void initialize();

        MeshData& ensureMesh(PrimitiveType type);
        MeshData& ensureMesh(const std::string& meshPath);
        MeshData& ensureMesh(const cardillo::MeshAsset& meshAsset);

    private:
        filament::Engine& engine_;
        uint32_t nextUniqueKey_ = 1;

        std::unordered_map<PrimitiveType, std::unique_ptr<MeshData>> primitiveMeshes_;
        std::unordered_map<std::string, std::unique_ptr<MeshData>> loadedMeshes_;
        std::unordered_map<const cardillo::MeshAsset*, std::unique_ptr<MeshData>> meshAssetMeshes_;
    };
};