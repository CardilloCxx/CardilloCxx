#include "geometry/Geometry.h"

#ifdef NORMAL
#undef NORMAL
#endif

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace frontend::geometry {

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
    if (length > 1e-6f) tangent /= length;
    float3 bitangent = cross(normal, tangent);
    filament::math::mat3f frame{tangent, bitangent, normal};
    return filament::math::mat3f::packTangentFrame(frame);
}

struct ObjVertexKey {
    int32_t pos = -1;
    int32_t nrm = -1;
    bool operator==(ObjVertexKey const& rhs) const noexcept {
        return pos == rhs.pos && nrm == rhs.nrm;
    }
};

struct ObjVertexKeyHash {
    size_t operator()(ObjVertexKey const& key) const noexcept {
        return (static_cast<size_t>(static_cast<uint32_t>(key.pos)) << 32) | static_cast<uint32_t>(key.nrm);
    }
};

MeshData loadObj(filament::Engine& engine, const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ObjLoader: failed to open " << path << "\n";
        return {};
    }

    std::vector<std::array<float,3>> positions;
    std::vector<std::array<float,3>> normals;
    std::vector<uint16_t> indices;
    std::unordered_map<ObjVertexKey, uint32_t, ObjVertexKeyHash> vertexMap;
    struct GVertex { float position[3]; float normal[3]; float tangent[4]; uint32_t color; };
    std::vector<GVertex> gverts;
    gverts.reserve(1024);

    auto fixIndex = [](int32_t index, size_t size) {
        return index < 0 ? static_cast<int32_t>(size) + index + 1 : index;
    };

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if (prefix == "v") {
            float x,y,z;
            iss >> x >> y >> z;
            positions.push_back({x,y,z});
        } else if (prefix == "vn") {
            float x,y,z;
            iss >> x >> y >> z;
            float len = std::sqrt(x*x + y*y + z*z);
            if (len > 1e-6f) {
                x /= len; y /= len; z /= len;
            }
            normals.push_back({x,y,z});
        } else if (prefix == "f") {
            struct FaceVertex { int32_t pos = -1; int32_t nrm = -1; };
            std::vector<FaceVertex> faceVerts;
            std::string token;
            bool foundMissingNormal = false;
            while (iss >> token) {
                FaceVertex fv;
                size_t firstSlash = token.find('/');
                if (firstSlash == std::string::npos) {
                    fv.pos = std::stoi(token);
                } else {
                    fv.pos = std::stoi(token.substr(0, firstSlash));
                    size_t secondSlash = token.find('/', firstSlash + 1);
                    if (secondSlash != std::string::npos && secondSlash + 1 < token.size()) {
                        fv.nrm = std::stoi(token.substr(secondSlash + 1));
                    }
                }
                fv.pos = fixIndex(fv.pos, positions.size()) - 1;
                if (fv.nrm != -1) {
                    fv.nrm = fixIndex(fv.nrm, normals.size()) - 1;
                } else if (!normals.empty()) {
                    foundMissingNormal = true;
                }
                faceVerts.push_back(fv);
            }
            if (faceVerts.size() >= 3) {
                for (size_t i = 1; i + 1 < faceVerts.size(); ++i) {
                    const FaceVertex tri[3] = { faceVerts[0], faceVerts[i], faceVerts[i + 1] };
                    for (int corner = 0; corner < 3; ++corner) {
                        const ObjVertexKey key{tri[corner].pos, tri[corner].nrm};
                        auto it = vertexMap.find(key);
                        if (it == vertexMap.end()) {
                            GVertex gv;
                            gv.position[0] = positions[key.pos][0];
                            gv.position[1] = positions[key.pos][1];
                            gv.position[2] = positions[key.pos][2];
                            if (key.nrm >= 0 && key.nrm < static_cast<int32_t>(normals.size())) {
                                gv.normal[0] = normals[key.nrm][0];
                                gv.normal[1] = normals[key.nrm][1];
                                gv.normal[2] = normals[key.nrm][2];
                            } else {
                                gv.normal[0] = 0.0f;
                                gv.normal[1] = 0.0f;
                                gv.normal[2] = 0.0f;
                            }
                            gv.tangent[0] = 0.0f;
                            gv.tangent[1] = 0.0f;
                            gv.tangent[2] = 1.0f;
                            gv.tangent[3] = 0.0f;
                            gv.color = 0xffffffffu;
                            const uint32_t newIndex = static_cast<uint32_t>(gverts.size());
                            vertexMap.emplace(key, newIndex);
                            gverts.push_back(gv);
                            it = vertexMap.find(key);
                        }
                        indices.push_back(static_cast<uint16_t>(it->second));
                    }
                }
                if (!normals.empty() && foundMissingNormal) {
                    for (auto& gv : gverts) {
                        gv.normal[0] = 0.0f;
                        gv.normal[1] = 0.0f;
                        gv.normal[2] = 0.0f;
                    }
                }
            }
        }
    }

    if (gverts.empty() || indices.empty()) return {};

    const bool haveNormals = !normals.empty();
    if (!haveNormals) {
        std::vector<filament::math::float3> accumulated(gverts.size(), filament::math::float3{0.0f, 0.0f, 0.0f});
        for (size_t i = 0; i + 2 < indices.size(); i += 3) {
            const uint16_t i0 = indices[i];
            const uint16_t i1 = indices[i + 1];
            const uint16_t i2 = indices[i + 2];
            filament::math::float3 p0{gverts[i0].position[0], gverts[i0].position[1], gverts[i0].position[2]};
            filament::math::float3 p1{gverts[i1].position[0], gverts[i1].position[1], gverts[i1].position[2]};
            filament::math::float3 p2{gverts[i2].position[0], gverts[i2].position[1], gverts[i2].position[2]};
            filament::math::float3 edge1 = p1 - p0;
            filament::math::float3 edge2 = p2 - p0;
            filament::math::float3 faceNormal = cross(edge2, edge1);
            accumulated[i0] += faceNormal;
            accumulated[i1] += faceNormal;
            accumulated[i2] += faceNormal;
        }
        for (size_t vi = 0; vi < gverts.size(); ++vi) {
            filament::math::float3 normal = accumulated[vi];
            float len = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            if (len > 1e-6f) {
                normal /= len;
            } else {
                normal = filament::math::float3{0.0f, 1.0f, 0.0f};
            }
            gverts[vi].normal[0] = normal.x;
            gverts[vi].normal[1] = normal.y;
            gverts[vi].normal[2] = normal.z;
        }
    }

    for (auto& gv : gverts) {
        filament::math::float3 normal{gv.normal[0], gv.normal[1], gv.normal[2]};
        const filament::math::quatf tangentFrame = computeTangentFrame(normal);
        gv.tangent[0] = tangentFrame.x;
        gv.tangent[1] = tangentFrame.y;
        gv.tangent[2] = tangentFrame.z;
        gv.tangent[3] = tangentFrame.w;
    }

    MeshData mb;
    auto vb = filament::VertexBuffer::Builder()
        .vertexCount(static_cast<uint32_t>(gverts.size()))
        .bufferCount(1)
        .attribute(filament::VertexAttribute::POSITION, 0, filament::VertexBuffer::AttributeType::FLOAT3, offsetof(GVertex, position), sizeof(GVertex))
        .attribute(filament::VertexAttribute::TANGENTS, 0, filament::VertexBuffer::AttributeType::FLOAT4, offsetof(GVertex, tangent), sizeof(GVertex))
        .attribute(filament::VertexAttribute::COLOR, 0, filament::VertexBuffer::AttributeType::UBYTE4, offsetof(GVertex, color), sizeof(GVertex))
        .normalized(filament::VertexAttribute::COLOR)
        .build(engine);
    size_t vbSize = sizeof(GVertex) * gverts.size();
    uint8_t* vbCopy = new uint8_t[vbSize];
    std::memcpy(vbCopy, gverts.data(), vbSize);
    vb->setBufferAt(engine, 0, filament::VertexBuffer::BufferDescriptor(vbCopy, vbSize, [](void* buffer, size_t, void*){ delete[] reinterpret_cast<uint8_t*>(buffer); }, nullptr));

    auto ib = filament::IndexBuffer::Builder()
        .indexCount(static_cast<uint32_t>(indices.size()))
        .bufferType(filament::IndexBuffer::IndexType::USHORT)
        .build(engine);
    size_t ibSize = sizeof(uint16_t) * indices.size();
    uint8_t* ibCopy = new uint8_t[ibSize];
    std::memcpy(ibCopy, indices.data(), ibSize);
    ib->setBuffer(engine, filament::IndexBuffer::BufferDescriptor(ibCopy, ibSize, [](void* buffer, size_t, void*){ delete[] reinterpret_cast<uint8_t*>(buffer); }, nullptr));

    mb.vb = vb;
    mb.ib = ib;
    mb.indexCount = static_cast<uint32_t>(indices.size());
    mb.primitiveType = filament::RenderableManager::PrimitiveType::TRIANGLES;

    filament::math::float3 mn{gverts[0].position[0], gverts[0].position[1], gverts[0].position[2]};
    filament::math::float3 mx = mn;
    for (size_t i = 1; i < gverts.size(); ++i) {
        mn.x = std::min(mn.x, gverts[i].position[0]); mn.y = std::min(mn.y, gverts[i].position[1]); mn.z = std::min(mn.z, gverts[i].position[2]);
        mx.x = std::max(mx.x, gverts[i].position[0]); mx.y = std::max(mx.y, gverts[i].position[1]); mx.z = std::max(mx.z, gverts[i].position[2]);
    }
    mb.minPos = mn;
    mb.maxPos = mx;
    return mb;
}

} // namespace filament_example::geometry
