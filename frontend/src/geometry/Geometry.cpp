#include "geometry/Geometry.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <cstddef>
#include <math/vec3.h>
#include <math/mat3.h>
#include <math/quat.h>

#ifdef NORMAL
#undef NORMAL
#endif

namespace {
    static void bufferFreeCallback(void* buffer, size_t /*size*/, void* /*user*/) {
        delete[] reinterpret_cast<uint8_t*>(buffer);
    }
}

namespace frontend::geometry {

struct Vertex {
    float position[3];
    float normal[3];
    float tangent[4];
    float uv[2];
    uint32_t color;
};

static uint32_t packColor(uint8_t r, uint8_t g, uint8_t b, uint8_t a = 0xff) {
    return (static_cast<uint32_t>(a) << 24) |
           (static_cast<uint32_t>(b) << 16) |
           (static_cast<uint32_t>(g) << 8)  |
           static_cast<uint32_t>(r);
}

static MeshData buildBuffersFromData(filament::Engine& engine, const std::vector<Vertex>& vertices, const std::vector<uint16_t>& indices, const std::vector<filament::math::float2>& /*uvs*/) {
    MeshData mb;
    if (vertices.empty() || indices.empty()) return mb;

    mb.vb = filament::VertexBuffer::Builder()
        .vertexCount(static_cast<uint32_t>(vertices.size()))
        .bufferCount(1)
        .attribute(filament::VertexAttribute::POSITION, 0, filament::VertexBuffer::AttributeType::FLOAT3, offsetof(Vertex, position), sizeof(Vertex))
        .attribute(filament::VertexAttribute::TANGENTS, 0, filament::VertexBuffer::AttributeType::FLOAT4, offsetof(Vertex, tangent), sizeof(Vertex))
        .attribute(filament::VertexAttribute::UV0, 0, filament::VertexBuffer::AttributeType::FLOAT2, offsetof(Vertex, uv), sizeof(Vertex))
        .build(engine);

    // Allocate a persistent copy for the VertexBuffer and provide a deleter callback
    {
        const size_t vbSize = sizeof(Vertex) * vertices.size();
        uint8_t* vbCopy = new uint8_t[vbSize];
        std::memcpy(vbCopy, vertices.data(), vbSize);
        mb.vb->setBufferAt(engine, 0, filament::VertexBuffer::BufferDescriptor(vbCopy, vbSize, bufferFreeCallback, nullptr));
    }

    mb.ib = filament::IndexBuffer::Builder()
        .indexCount(static_cast<uint32_t>(indices.size()))
        .bufferType(filament::IndexBuffer::IndexType::USHORT)
        .build(engine);

    // Allocate a persistent copy for the IndexBuffer and provide a deleter callback
    {
        const size_t ibSize = sizeof(uint16_t) * indices.size();
        uint8_t* ibCopy = new uint8_t[ibSize];
        std::memcpy(ibCopy, indices.data(), ibSize);
        mb.ib->setBuffer(engine, filament::IndexBuffer::BufferDescriptor(ibCopy, ibSize, bufferFreeCallback, nullptr));
    }

    mb.indexCount = static_cast<uint32_t>(indices.size());
    mb.primitiveType = filament::RenderableManager::PrimitiveType::TRIANGLES;

    // Compute min/max bounds
    filament::math::float3 mn{vertices[0].position[0], vertices[0].position[1], vertices[0].position[2]};
    filament::math::float3 mx = mn;
    for (size_t i = 1; i < vertices.size(); ++i) {
        mn.x = std::min(mn.x, vertices[i].position[0]);
        mn.y = std::min(mn.y, vertices[i].position[1]);
        mn.z = std::min(mn.z, vertices[i].position[2]);
        mx.x = std::max(mx.x, vertices[i].position[0]);
        mx.y = std::max(mx.y, vertices[i].position[1]);
        mx.z = std::max(mx.z, vertices[i].position[2]);
    }
    mb.minPos = mn;
    mb.maxPos = mx;
    mb.bounds = filament::Box{};
    mb.bounds.set(mn, mx);

    return mb;
}

static filament::math::quatf computeTangentFrame(const filament::math::float3& normal) {
    using filament::math::float3;
    float3 up = std::fabs(normal.x) > 0.9f ? float3{0.0f, 1.0f, 0.0f} : float3{1.0f, 0.0f, 0.0f};
    float3 tangent = cross(up, normal);
    float length = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z);
    if (length < 1e-6f) {
        up = float3{0.0f, 0.0f, 1.0f};
        tangent = cross(up, normal);
        length = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z);
    }
    if (length > 1e-6f) {
        tangent /= length;
    }
    float3 bitangent = cross(normal, tangent);
    filament::math::mat3f frame{tangent, bitangent, normal};
    return filament::math::mat3f::packTangentFrame(frame);
}

static bool computeTriangleTangentBasis(const filament::math::float3& p0,
        const filament::math::float3& p1,
        const filament::math::float3& p2,
        const filament::math::float2& uv0,
        const filament::math::float2& uv1,
        const filament::math::float2& uv2,
        filament::math::float3& tangent,
        filament::math::float3& bitangent) {
    using filament::math::float2;
    const filament::math::float3 edge1 = p1 - p0;
    const filament::math::float3 edge2 = p2 - p0;
    const float2 deltaUV1 = uv1 - uv0;
    const float2 deltaUV2 = uv2 - uv0;
    const float denom = deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y;
    if (std::fabs(denom) < 1e-8f) {
        return false;
    }
    const float invDet = 1.0f / denom;
    tangent = edge1 * deltaUV2.y - edge2 * deltaUV1.y;
    bitangent = edge2 * deltaUV1.x - edge1 * deltaUV2.x;
    tangent *= invDet;
    bitangent *= invDet;
    return true;
}

static filament::math::quatf computeTangentFrameForTriangle(const filament::math::float3& p0,
        const filament::math::float3& p1,
        const filament::math::float3& p2,
        const filament::math::float2& uv0,
        const filament::math::float2& uv1,
        const filament::math::float2& uv2,
        const filament::math::float3& normal) {
    using filament::math::float3;
    filament::math::float3 tangent;
    filament::math::float3 bitangent;
    if (computeTriangleTangentBasis(p0, p1, p2, uv0, uv1, uv2, tangent, bitangent)) {
        const float tlen = std::sqrt(dot(tangent, tangent));
        if (tlen > 1e-6f) {
            tangent /= tlen;
            float3 bitan = cross(normal, tangent);
            const float handedness = (dot(bitan, bitangent) < 0.0f) ? -1.0f : 1.0f;
            bitan *= handedness;
            filament::math::mat3f frame{tangent, bitan, normal};
            return filament::math::mat3f::packTangentFrame(frame);
        }
    }
    return computeTangentFrame(normal);
}

static MeshData buildFlatFromData(filament::Engine& engine, const std::vector<Vertex>& vertices, const std::vector<uint16_t>& indices, const std::vector<filament::math::float2>& uvs) {
    using filament::math::float3;

    std::vector<Vertex> flatVertices;
    std::vector<uint16_t> flatIndices;
    flatVertices.reserve(indices.size());
    flatIndices.reserve(indices.size());

    for (size_t i = 0; i < indices.size(); i += 3) {
        const uint16_t i0 = indices[i + 0];
        const uint16_t i1 = indices[i + 1];
        const uint16_t i2 = indices[i + 2];

        const float3 p0{vertices[i0].position[0], vertices[i0].position[1], vertices[i0].position[2]};
        const float3 p1{vertices[i1].position[0], vertices[i1].position[1], vertices[i1].position[2]};
        const float3 p2{vertices[i2].position[0], vertices[i2].position[1], vertices[i2].position[2]};

        float3 normal = cross(p1 - p0, p2 - p0);
        const float len = length(normal);
        normal = (len > 1e-6f) ? (normal / len) : float3{0.0f, 0.0f, 1.0f};

        const filament::math::float2 uv0 = (i0 < uvs.size()) ? uvs[i0] : filament::math::float2{vertices[i0].uv[0], vertices[i0].uv[1]};
        const filament::math::float2 uv1 = (i1 < uvs.size()) ? uvs[i1] : filament::math::float2{vertices[i1].uv[0], vertices[i1].uv[1]};
        const filament::math::float2 uv2 = (i2 < uvs.size()) ? uvs[i2] : filament::math::float2{vertices[i2].uv[0], vertices[i2].uv[1]};
        const filament::math::quatf tangentFrame = computeTangentFrameForTriangle(p0, p1, p2, uv0, uv1, uv2, normal);

        for (int corner = 0; corner < 3; ++corner) {
            const uint16_t origIdx = indices[i + corner];
            Vertex dst = vertices[origIdx];
            dst.normal[0] = normal.x;
            dst.normal[1] = normal.y;
            dst.normal[2] = normal.z;
            dst.tangent[0] = tangentFrame.x;
            dst.tangent[1] = tangentFrame.y;
            dst.tangent[2] = tangentFrame.z;
            dst.tangent[3] = tangentFrame.w;
            if (origIdx < uvs.size()) {
                dst.uv[0] = uvs[origIdx].x;
                dst.uv[1] = uvs[origIdx].y;
            }
            flatVertices.push_back(dst);
            flatIndices.push_back(static_cast<uint16_t>(flatVertices.size() - 1));
        }
    }

    return buildBuffersFromData(engine, flatVertices, flatIndices, uvs);
}

static MeshData buildSmoothFromData(filament::Engine& engine, std::vector<Vertex> vertices, const std::vector<uint16_t>& indices) {
    using filament::math::float3;

    std::vector<float3> normals(vertices.size(), float3{0.0f});
    std::vector<float3> tangents(vertices.size(), float3{0.0f});
    std::vector<float3> bitangents(vertices.size(), float3{0.0f});
    std::vector<filament::math::float2> uvs(vertices.size(), filament::math::float2{0.0f, 0.0f});

    for (size_t i = 0; i < indices.size(); i += 3) {
        const uint16_t i0 = indices[i + 0];
        const uint16_t i1 = indices[i + 1];
        const uint16_t i2 = indices[i + 2];

        const float3 p0{vertices[i0].position[0], vertices[i0].position[1], vertices[i0].position[2]};
        const float3 p1{vertices[i1].position[0], vertices[i1].position[1], vertices[i1].position[2]};
        const float3 p2{vertices[i2].position[0], vertices[i2].position[1], vertices[i2].position[2]};

        const float3 faceNormal = cross(p1 - p0, p2 - p0);

        normals[i0] += faceNormal;
        normals[i1] += faceNormal;
        normals[i2] += faceNormal;

        const filament::math::float2 uv0{vertices[i0].uv[0], vertices[i0].uv[1]};
        const filament::math::float2 uv1{vertices[i1].uv[0], vertices[i1].uv[1]};
        const filament::math::float2 uv2{vertices[i2].uv[0], vertices[i2].uv[1]};

        filament::math::float3 triTangent;
        filament::math::float3 triBitangent;
        if (computeTriangleTangentBasis(p0, p1, p2, uv0, uv1, uv2, triTangent, triBitangent)) {
            tangents[i0] += triTangent;
            tangents[i1] += triTangent;
            tangents[i2] += triTangent;

            bitangents[i0] += triBitangent;
            bitangents[i1] += triBitangent;
            bitangents[i2] += triBitangent;
        }
    }

    for (size_t vi = 0; vi < vertices.size(); ++vi) {
        float3 normal = normals[vi];
        const float len = length(normal);
        normal = (len > 1e-6f) ? (normal / len) : float3{0.0f, 1.0f, 0.0f};

        float3 tangent = tangents[vi];
        const float tanLen = std::sqrt(dot(tangent, tangent));
        if (tanLen > 1e-6f) {
            tangent /= tanLen;
            tangent -= normal * dot(normal, tangent);
            const float correctedLen = std::sqrt(dot(tangent, tangent));
            if (correctedLen > 1e-6f) {
                tangent /= correctedLen;
            } else {
                float3 up = std::fabs(normal.x) > 0.9f ? float3{0.0f, 1.0f, 0.0f} : float3{1.0f, 0.0f, 0.0f};
                tangent = normalize(cross(up, normal));
            }
        } else {
            float3 up = std::fabs(normal.x) > 0.9f ? float3{0.0f, 1.0f, 0.0f} : float3{1.0f, 0.0f, 0.0f};
            tangent = normalize(cross(up, normal));
        }

        float3 bitangent = cross(normal, tangent);
        if (length(bitangents[vi]) > 1e-6f) {
            const float handedness = (dot(bitangent, bitangents[vi]) < 0.0f) ? -1.0f : 1.0f;
            bitangent *= handedness;
        }

        filament::math::mat3f frame{tangent, bitangent, normal};
        const filament::math::quatf tangentFrame = filament::math::mat3f::packTangentFrame(frame);

        vertices[vi].normal[0] = normal.x;
        vertices[vi].normal[1] = normal.y;
        vertices[vi].normal[2] = normal.z;
        vertices[vi].tangent[0] = tangentFrame.x;
        vertices[vi].tangent[1] = tangentFrame.y;
        vertices[vi].tangent[2] = tangentFrame.z;
        vertices[vi].tangent[3] = tangentFrame.w;

        uvs[vi] = filament::math::float2{vertices[vi].uv[0], vertices[vi].uv[1]};
    }

    return buildBuffersFromData(engine, vertices, indices, uvs);
}

MeshData buildMeshFromTriangles(filament::Engine& engine,
        const std::vector<filament::math::float3>& positions,
        const std::vector<uint16_t>& indices) {
    std::vector<Vertex> vertices;
    vertices.reserve(positions.size());
    for (const auto& position : positions) {
        Vertex v{};
        v.position[0] = position.x;
        v.position[1] = position.y;
        v.position[2] = position.z;
        v.color = packColor(0xff, 0xff, 0xff);
        vertices.push_back(v);
    }

    return buildSmoothFromData(engine, vertices, indices);
}

MeshData buildCube(filament::Engine& engine) {
    std::vector<Vertex> vertices;
    std::vector<uint16_t> indices;
    vertices.reserve(24);
    indices.reserve(36);

    const uint32_t white = 0xffffffffu;
    const float s = 0.5f;

    // Define 6 faces with positions (4 per face) and standard quad UVs
    struct FaceDef {
        float pos[4][3];
        float uv[4][2];
    };

    const FaceDef faces[6] = {
        // Back face (z = -s)
        {{{ s, -s, -s}, {-s, -s, -s}, {-s,  s, -s}, { s,  s, -s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}},
        // Front face (z = +s)
        {{{-s, -s,  s}, { s, -s,  s}, { s,  s,  s}, {-s,  s,  s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}},
        // Bottom face (y = -s)
        {{{-s, -s, -s}, { s, -s, -s}, { s, -s,  s}, {-s, -s,  s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}},
        // Top face (y = +s)
        {{{-s,  s,  s}, { s,  s,  s}, { s,  s, -s}, {-s,  s, -s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}},
        // Left face (x = -s)
        {{{-s, -s, -s}, {-s, -s,  s}, {-s,  s,  s}, {-s,  s, -s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}},
        // Right face (x = +s)
        {{{ s, -s,  s}, { s, -s, -s}, { s,  s, -s}, { s,  s,  s}}, {{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}}
    };

    for (int f = 0; f < 6; ++f) {
        const uint16_t baseIdx = static_cast<uint16_t>(vertices.size());

        for (int v = 0; v < 4; ++v) {
            Vertex vert{};
            vert.position[0] = faces[f].pos[v][0];
            vert.position[1] = faces[f].pos[v][1];
            vert.position[2] = faces[f].pos[v][2];
            vert.uv[0]       = faces[f].uv[v][0];
            vert.uv[1]       = faces[f].uv[v][1];
            vert.color       = white;
            vertices.push_back(vert);
        }

        // Two triangles per face quad (CCW winding)
        indices.push_back(baseIdx + 0);
        indices.push_back(baseIdx + 1);
        indices.push_back(baseIdx + 2);

        indices.push_back(baseIdx + 0);
        indices.push_back(baseIdx + 2);
        indices.push_back(baseIdx + 3);
    }

    std::vector<filament::math::float2> uvs; // UVs are pre-populated on vertices
    return buildFlatFromData(engine, vertices, indices, uvs);
}

MeshData buildSphere(filament::Engine& engine, int lat, int lon) {
    std::vector<Vertex> vertices;
    std::vector<uint16_t> indices;
    const uint32_t white = packColor(0xff, 0xff, 0xff);

    for (int i = 0; i <= lat; ++i) {
        const float v = static_cast<float>(i) / lat;
        const float theta = v * 3.14159265358979323846f;
        const float sinTheta = std::sin(theta);
        const float cosTheta = std::cos(theta);

        for (int j = 0; j <= lon; ++j) {
            const float u = static_cast<float>(j) / lon;
            const float phi = u * 2.0f * 3.14159265358979323846f;
            const float sinPhi = std::sin(phi);
            const float cosPhi = std::cos(phi);

            const float x = sinTheta * cosPhi * 0.5f;
            const float y = cosTheta * 0.5f;
            const float z = sinTheta * sinPhi * 0.5f;

            // 1. Exact Outward Normal
            filament::math::float3 normal = normalize(filament::math::float3{x, y, z});

            // 2. Exact Tangent Along Latitude
            filament::math::float3 tangent{-sinPhi, 0.0f, cosPhi};
            if (std::abs(sinTheta) < 1e-5f) {
                tangent = filament::math::float3{1.0f, 0.0f, 0.0f}; // Pole fallback
            } else {
                tangent = normalize(tangent);
            }

            // 3. Orthonormal Basis & Tangent Frame Quaternion
            filament::math::float3 bitangent = cross(normal, tangent);
            filament::math::mat3f frame{tangent, bitangent, normal};
            filament::math::quatf q = filament::math::mat3f::packTangentFrame(frame);

            Vertex vert{};
            vert.position[0] = x;
            vert.position[1] = y;
            vert.position[2] = z;
            vert.normal[0] = normal.x;
            vert.normal[1] = normal.y;
            vert.normal[2] = normal.z;
            vert.tangent[0] = q.x;
            vert.tangent[1] = q.y;
            vert.tangent[2] = q.z;
            vert.tangent[3] = q.w;
            vert.uv[0] = u;
            vert.uv[1] = v;
            vert.color = white;

            vertices.push_back(vert);
        }
    }

    for (int i = 0; i < lat; ++i) {
        for (int j = 0; j < lon; ++j) {
            const uint16_t a = static_cast<uint16_t>(i * (lon + 1) + j);
            const uint16_t b = static_cast<uint16_t>(a + lon + 1);

            indices.push_back(a);
            indices.push_back(a + 1);
            indices.push_back(b);

            indices.push_back(a + 1);
            indices.push_back(b + 1);
            indices.push_back(b);
        }
    }

    std::vector<filament::math::float2> emptyUvs;
    return buildBuffersFromData(engine, vertices, indices, emptyUvs);
}

MeshData buildCylinder(filament::Engine& engine, int seg) {
    std::vector<Vertex> vertices;
    std::vector<uint16_t> indices;
    const uint32_t white = packColor(0xff, 0xff, 0xff);
    
    constexpr float radius = 0.5f;
    constexpr float height = 1.0f;
    const float halfHeight = height * 0.5f;

    for (int i = 0; i < seg; ++i) {
        const float u = static_cast<float>(i) / seg;
        const float a = u * 2.0f * 3.14159265f;
        const float x = std::cos(a) * radius;
        const float z = std::sin(a) * radius;

        Vertex vt{};
        vt.position[0] = x; vt.position[1] = halfHeight; vt.position[2] = z;
        vt.uv[0] = u; vt.uv[1] = 1.0f; vt.color = white;
        vertices.push_back(vt);

        Vertex vb{};
        vb.position[0] = x; vb.position[1] = -halfHeight; vb.position[2] = z;
        vb.uv[0] = u; vb.uv[1] = 0.0f; vb.color = white;
        vertices.push_back(vb);
    }

    for (int i = 0; i < seg; ++i) {
        const uint16_t topA = static_cast<uint16_t>(i * 2);
        const uint16_t botA = static_cast<uint16_t>(i * 2 + 1);
        const uint16_t topB = static_cast<uint16_t>(((i + 1) % seg) * 2);
        const uint16_t botB = static_cast<uint16_t>(((i + 1) % seg) * 2 + 1);

        indices.push_back(topA); indices.push_back(topB); indices.push_back(botA);
        indices.push_back(topB); indices.push_back(botB); indices.push_back(botA);
    }

    const uint16_t topCenter = static_cast<uint16_t>(vertices.size());
    Vertex topVertex{};
    topVertex.position[0] = 0.0f; topVertex.position[1] = halfHeight; topVertex.position[2] = 0.0f;
    topVertex.uv[0] = 0.5f; topVertex.uv[1] = 0.5f;
    topVertex.color = white; vertices.push_back(topVertex);

    const uint16_t bottomCenter = static_cast<uint16_t>(vertices.size());
    Vertex bottomVertex{};
    bottomVertex.position[0] = 0.0f; bottomVertex.position[1] = -halfHeight; bottomVertex.position[2] = 0.0f;
    bottomVertex.uv[0] = 0.5f; bottomVertex.uv[1] = 0.5f;
    bottomVertex.color = white; vertices.push_back(bottomVertex);

    for (int i = 0; i < seg; ++i) {
        const uint16_t topA = static_cast<uint16_t>(i * 2);
        const uint16_t topB = static_cast<uint16_t>(((i + 1) % seg) * 2);
        
        // Top cap
        indices.push_back(topCenter); indices.push_back(topB); indices.push_back(topA);

        const uint16_t botA = static_cast<uint16_t>(i * 2 + 1);
        const uint16_t botB = static_cast<uint16_t>(((i + 1) % seg) * 2 + 1);
        
        // Bottom cap
        indices.push_back(bottomCenter); indices.push_back(botA); indices.push_back(botB);
    }

    return buildSmoothFromData(engine, vertices, indices);
}

MeshData buildCone(filament::Engine& engine, int seg) {
    std::vector<Vertex> vertices;
    std::vector<uint16_t> indices;
    const uint32_t white = packColor(0xff, 0xff, 0xff);
    const float radius = 0.5f;
    const float height = 1.0f;
    for (int i = 0; i <= seg; ++i) {
        const float u = static_cast<float>(i) / seg;
        const float a = u * 2.0f * 3.14159265f;
        const float x = std::cos(a) * radius;
        const float z = std::sin(a) * radius;
        Vertex vb{};
        vb.position[0] = x;
        vb.position[1] = -height * 0.5f;
        vb.position[2] = z;
        vb.uv[0] = u;
        vb.uv[1] = 0.0f;
        vb.color = white;
        vertices.push_back(vb);
    }
    // apex
    Vertex apex{};
    apex.position[0] = 0.0f;
    apex.position[1] = height * 0.5f;
    apex.position[2] = 0.0f;
    apex.uv[0] = 0.5f;
    apex.uv[1] = 1.0f;
    apex.color = white;
    vertices.push_back(apex);

    const uint16_t apexIdx = static_cast<uint16_t>(vertices.size() - 1);
    const int base = 0;
    for (int i = 0; i < seg; ++i) {
        const uint16_t a = static_cast<uint16_t>(base + i);
        const uint16_t b = static_cast<uint16_t>(base + ((i + 1) % seg));
        indices.push_back(apexIdx); indices.push_back(b); indices.push_back(a);
    }
    return buildSmoothFromData(engine, vertices, indices);
}

MeshData buildCapsule(filament::Engine& engine) {
    std::vector<Vertex> vertices;
    std::vector<uint16_t> indices;
    
    constexpr int segments = 24;
    constexpr int capRings = 6;
    constexpr float radius = 0.5f;
    constexpr float halfLength = 0.5f;

    for (int ring = 0; ring <= capRings; ++ring) {
        float theta = -3.14159265f / 2.0f + (static_cast<float>(ring) / capRings) * (3.14159265f / 2.0f);
        float r = radius * std::cos(theta);
        float z = -halfLength + radius * std::sin(theta);

        for (int i = 0; i < segments; ++i) {
            float u = static_cast<float>(i) / segments;
            float angle = u * 2.0f * 3.14159265f;
            Vertex v{};
            v.position[0] = std::cos(angle) * r;
            v.position[1] = z;
            v.position[2] = std::sin(angle) * r;
            v.uv[0] = u;
            v.uv[1] = (z + 1.0f) * 0.5f;
            vertices.push_back(v);
        }
    }

    for (int ring = 0; ring <= capRings; ++ring) {
        float theta = (static_cast<float>(ring) / capRings) * (3.14159265f / 2.0f);
        float r = radius * std::cos(theta);
        float z = halfLength + radius * std::sin(theta);

        for (int i = 0; i < segments; ++i) {
            float u = static_cast<float>(i) / segments;
            float angle = u * 2.0f * 3.14159265f;
            Vertex v{};
            v.position[0] = std::cos(angle) * r;
            v.position[1] = z;
            v.position[2] = std::sin(angle) * r;
            v.uv[0] = u;
            v.uv[1] = (z + 1.0f) * 0.5f;
            vertices.push_back(v);
        }
    }

    int totalRings = (capRings + 1) * 2;
    for (int ring = 0; ring + 1 < totalRings; ++ring) {
        uint16_t ring0Base = ring * segments;
        uint16_t ring1Base = (ring + 1) * segments;

        for (int i = 0; i < segments; ++i) {
            uint16_t next = (i + 1) % segments;
            
            uint16_t a = ring0Base + i;
            uint16_t b = ring0Base + next;
            uint16_t c = ring1Base + i;
            uint16_t d = ring1Base + next;

            indices.push_back(a);
            indices.push_back(b);
            indices.push_back(c);

            indices.push_back(c);
            indices.push_back(b);
            indices.push_back(d);
        }
    }

    return buildSmoothFromData(engine, vertices, indices);
}

MeshData buildTestTriangle(filament::Engine& engine) {
    struct TV { float pos[3]; uint32_t color; };
    std::vector<TV> verts(3);
    std::vector<uint16_t> idx;
    verts[0].pos[0] = -0.5f; verts[0].pos[1] = -0.5f; verts[0].pos[2] = 0.0f; verts[0].color = packColor(0xff,0xff,0xff);
    verts[1].pos[0] =  0.5f; verts[1].pos[1] = -0.5f; verts[1].pos[2] = 0.0f; verts[1].color = packColor(0xff,0xff,0xff);
    verts[2].pos[0] =  0.0f; verts[2].pos[1] =  0.5f; verts[2].pos[2] = 0.0f; verts[2].color = packColor(0xff,0xff,0xff);
    idx.push_back(0); idx.push_back(1); idx.push_back(2);

    std::vector<Vertex> outVerts;
    outVerts.reserve(3);
    for (size_t i = 0; i < verts.size(); ++i) {
        Vertex vv{};
        vv.position[0] = verts[i].pos[0];
        vv.position[1] = verts[i].pos[1];
        vv.position[2] = verts[i].pos[2];
        vv.color = verts[i].color;
        vv.uv[0] = (i == 1) ? 1.0f : 0.0f;
        vv.uv[1] = (i == 2) ? 1.0f : 0.0f;
        outVerts.push_back(vv);
    }
    return buildSmoothFromData(engine, outVerts, idx);
}

} // namespace frontend::geometry