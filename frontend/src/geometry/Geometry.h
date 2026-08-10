#pragma once

#ifdef NORMAL
#undef NORMAL
#endif

#include <filament/Engine.h>
#include <filament/IndexBuffer.h>
#include <filament/VertexBuffer.h>
#include <filament/RenderableManager.h>
#include <vector>

struct MeshData {
    filament::VertexBuffer* vb = nullptr;
    filament::IndexBuffer* ib = nullptr;
    uint32_t indexCount = 0;
    filament::RenderableManager::PrimitiveType primitiveType;
    filament::math::float3 minPos;
    filament::math::float3 maxPos;
    filament::Box bounds;

    uint32_t uniqueKey = 0;
};

namespace frontend::geometry {

MeshData buildCube(filament::Engine& engine);
MeshData buildSphere(filament::Engine& engine, int lat = 20, int lon = 24);
MeshData buildCylinder(filament::Engine& engine, int seg = 20);
MeshData buildCone(filament::Engine& engine, int seg = 20);
MeshData buildCapsule(filament::Engine& engine);
// Simple test triangle to ensure pipeline renders a known-good primitive.
MeshData buildTestTriangle(filament::Engine& engine);

// Load a simple wavefront OBJ (positions + faces). Returns empty MeshData on failure.
MeshData loadObj(filament::Engine& engine, const std::string& path);

// Build a triangle mesh directly from explicit vertices and indices.
MeshData buildMeshFromTriangles(filament::Engine& engine,
    const std::vector<filament::math::float3>& positions,
    const std::vector<uint16_t>& indices);

} // namespace frontend::geometry
