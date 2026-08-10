#include "geometry/Assets.h"
#include "physics/CardilloPhysicsBridge.h"

namespace frontend {
    Assets::Assets(filament::Engine& engine) : engine_(engine) {}

    Assets::~Assets() {
        for (auto& [type, mesh] : primitiveMeshes_) {
            if (mesh->vb) engine_.destroy(mesh->vb);
            if (mesh->ib) engine_.destroy(mesh->ib);
        }
        for (auto& [path, mesh] : loadedMeshes_) {
            if (mesh->vb) engine_.destroy(mesh->vb);
            if (mesh->ib) engine_.destroy(mesh->ib);
        }
        for (auto& [asset, mesh] : meshAssetMeshes_) {
            if (mesh->vb) engine_.destroy(mesh->vb);
            if (mesh->ib) engine_.destroy(mesh->ib);
        }
    }

    void Assets::initialize() {
        ensureMesh(PrimitiveType::Cube);
        ensureMesh(PrimitiveType::Sphere);
        ensureMesh(PrimitiveType::Capsule);
        ensureMesh(PrimitiveType::Cylinder);
        ensureMesh(PrimitiveType::Cone);
    }

    MeshData& Assets::ensureMesh(PrimitiveType type) {
        auto it = primitiveMeshes_.find(type);
        if (it != primitiveMeshes_.end()) return *it->second;

        MeshData mesh;
        switch (type) {
            case PrimitiveType::Cube: mesh = geometry::buildCube(engine_); break;
            case PrimitiveType::Sphere: mesh = geometry::buildSphere(engine_); break;
            case PrimitiveType::Capsule: mesh = geometry::buildCapsule(engine_); break;
            case PrimitiveType::Cylinder: mesh = geometry::buildCylinder(engine_); break;
            case PrimitiveType::Cone: mesh = geometry::buildCone(engine_); break;
            default: throw std::runtime_error("Unsupported primitive type");
        }

        mesh.uniqueKey = nextUniqueKey_++;
        auto inserted = primitiveMeshes_.emplace(type, std::make_unique<MeshData>(std::move(mesh)));
        return *inserted.first->second;
    }

    MeshData& Assets::ensureMesh(const std::string& meshPath) {
        auto it = loadedMeshes_.find(meshPath);
        if (it != loadedMeshes_.end()) return *it->second;

        auto mb = geometry::loadObj(engine_, meshPath);
        if (!mb.vb || !mb.ib) {
            throw std::runtime_error("Failed to load mesh: " + meshPath);
        }

        mb.uniqueKey = nextUniqueKey_++;
        auto inserted = loadedMeshes_.emplace(meshPath, std::make_unique<MeshData>(std::move(mb)));
        return *inserted.first->second;
    }

    MeshData& Assets::ensureMesh(const cardillo::MeshAsset& meshAsset) {
        auto it = meshAssetMeshes_.find(&meshAsset);
        if (it != meshAssetMeshes_.end()) return *it->second;

        std::vector<filament::math::float3> vertices;
        std::vector<uint16_t> indices;
        if (meshAsset.bvh && meshAsset.bvh->vertices && meshAsset.bvh->tri_indices) {
            const auto& rawVertices = *meshAsset.bvh->vertices;
            const auto& rawTriangles = *meshAsset.bvh->tri_indices;
            vertices.reserve(rawVertices.size());
            indices.reserve(rawTriangles.size() * 3);
            for (const auto& rawVertex : rawVertices) {
                const cardillo::Vector3r vertexPos(rawVertex[0], rawVertex[1], rawVertex[2]);
                vertices.emplace_back(physics::CardilloPhysicsBridge::toFloat3(vertexPos));
            }
            for (const auto& triangle : rawTriangles) {
                indices.push_back(static_cast<uint16_t>(triangle[0]));
                indices.push_back(static_cast<uint16_t>(triangle[1]));
                indices.push_back(static_cast<uint16_t>(triangle[2]));
            }
        }

        MeshData mesh = geometry::buildMeshFromTriangles(engine_, vertices, indices);
        mesh.uniqueKey = nextUniqueKey_++;

        auto inserted = meshAssetMeshes_.emplace(&meshAsset, std::make_unique<MeshData>(std::move(mesh)));
        return *inserted.first->second;
    }
}