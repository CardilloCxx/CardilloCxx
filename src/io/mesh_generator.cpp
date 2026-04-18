#include "mesh_generator.hpp"
#include <cmath>
#include "../rigid_body/rigid_body.hpp"

namespace cardillo::io {

void MeshGenerator::generateUVSphere(int latSegments, int lonSegments, std::vector<Vector3r>& outVertices, std::vector<Eigen::Vector3i>& outTriangles) {
    if (latSegments < 3) latSegments = 3;
    if (lonSegments < 3) lonSegments = 3;
    outVertices.clear();
    outTriangles.clear();
    // Top pole
    outVertices.emplace_back(0, 0, 1);
    // Rings
    for (int i = 1; i < latSegments; ++i) {
        const real_t phi = (real_t)M_PI * (real_t)i / (real_t)latSegments;  // (0,pi)
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
    outVertices.emplace_back(0, 0, -1);

    // Helper to index ring vertices
    auto ringIndex = [&](int ring, int seg) -> int { return 1 + ring * lonSegments + (seg % lonSegments); };

    // Top cap triangles
    for (int j = 0; j < lonSegments; ++j) {
        int a = 0;  // top pole
        int b = ringIndex(0, j);
        int c = ringIndex(0, j + 1);
        outTriangles.emplace_back(a, b, c);
    }
    // Middle quads split into triangles
    for (int i = 0; i < latSegments - 2; ++i) {
        for (int j = 0; j < lonSegments; ++j) {
            int v00 = ringIndex(i, j);
            int v01 = ringIndex(i, j + 1);
            int v10 = ringIndex(i + 1, j);
            int v11 = ringIndex(i + 1, j + 1);
            outTriangles.emplace_back(v00, v10, v11);
            outTriangles.emplace_back(v00, v11, v01);
        }
    }
    // Bottom cap triangles
    const int lastRing = latSegments - 2;
    for (int j = 0; j < lonSegments; ++j) {
        int a = ringIndex(lastRing, j);
        int b = bottomIndex;
        int c = ringIndex(lastRing, j + 1);
        outTriangles.emplace_back(a, b, c);
    }
}

void MeshGenerator::generateCapsuleMesh(int segmentsCircumference, int segmentsHemisphere, int segmentsCylinder, real_t radius, real_t halfLength, std::vector<Vector3r>& vertices,
                                        std::vector<Eigen::Vector3i>& triangles) {
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

void MeshGenerator::generateCylinderMesh(int segmentsCircumference, real_t radius, real_t halfLength, std::vector<Vector3r>& vertices, std::vector<Eigen::Vector3i>& triangles) {
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

void MeshGenerator::generateConeMesh(int segmentsCircumference, real_t radius, real_t height, std::vector<Vector3r>& vertices, std::vector<Eigen::Vector3i>& triangles) {
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

namespace {

enum class MeshKind { Unsupported, SoftBody, Plane, Cube, Capsule, Cylinder, Cone, MeshAsset, HeightField, Sphere };

MeshKind detectMeshKind(const entt::registry& reg, entt::entity e) {
    if (reg.any_of<cardillo::C_SoftBodySurface>(e)) return MeshKind::SoftBody;
    if (reg.any_of<cardillo::C_Plane>(e) && reg.any_of<cardillo::C_Position3>(e)) return MeshKind::Plane;
    if (reg.any_of<cardillo::C_Cube>(e)) return MeshKind::Cube;
    if (reg.any_of<cardillo::C_Capsule>(e)) return MeshKind::Capsule;
    if (reg.any_of<cardillo::C_Cylinder>(e)) return MeshKind::Cylinder;
    if (reg.any_of<cardillo::C_Cone>(e)) return MeshKind::Cone;
    if (reg.any_of<cardillo::C_Mesh>(e)) return MeshKind::MeshAsset;
    if (reg.any_of<cardillo::C_HeightField>(e)) return MeshKind::HeightField;
    if (reg.any_of<cardillo::C_RB_Sphere, cardillo::C_Radius>(e)) return MeshKind::Sphere;
    return MeshKind::Unsupported;
}

void fillCommonEntityData(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    out = cardillo::io::MeshGenerator::EntityMesh{};
    out.entity = e;
    out.entityId = static_cast<int>(entt::to_integral(e));
    out.isDynamic = reg.any_of<cardillo::C_PhysicsObject, cardillo::C_StaticTrajectory>(e);

    if (reg.any_of<cardillo::C_BeamElement>(e)) {
        const auto& be = reg.get<cardillo::C_BeamElement>(e);
        out.beamLengthRatio = (be.l0 > (real_t)0) ? (be.l / be.l0) : (real_t)1;
    }

    const auto state = cardillo::RigidBody::getState(reg, e);

    if (reg.any_of<cardillo::C_LinearVelocity3>(e)) {
        out.vlin = state.linearVelocity;
        out.hasKinematics = true;
    }
    if (reg.any_of<cardillo::C_AngularVelocity3>(e)) {
        out.omega = state.angularVelocity;
        out.hasKinematics = true;
    }
    if (reg.any_of<cardillo::C_LinearAcceleration3>(e)) {
        out.alin = reg.get<cardillo::C_LinearAcceleration3>(e).value;
    }
    if (reg.any_of<cardillo::C_AngularAcceleration3>(e)) {
        out.alpha = reg.get<cardillo::C_AngularAcceleration3>(e).value;
    }

    out.center = state.position;
    out.R = state.rotation;
}

bool buildSoftBodyMesh(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& surf = reg.get<cardillo::C_SoftBodySurface>(e);
    out.center.setZero();
    out.vertices.reserve(surf.nodes.size());
    out.perVertexVelocity.reserve(surf.nodes.size());
    out.perVertexAcceleration.reserve(surf.nodes.size());

    for (entt::entity node : surf.nodes) {
        if (reg.valid(node) && reg.any_of<cardillo::C_Position3>(node)) {
            out.vertices.push_back(reg.get<cardillo::C_Position3>(node).value);
        } else {
            out.vertices.emplace_back(Vector3r::Zero());
        }

        if (reg.valid(node) && reg.any_of<cardillo::C_LinearVelocity3>(node)) {
            out.perVertexVelocity.push_back(reg.get<cardillo::C_LinearVelocity3>(node).value);
        } else {
            out.perVertexVelocity.emplace_back(Vector3r::Zero());
        }

        if (reg.valid(node) && reg.any_of<cardillo::C_LinearAcceleration3>(node)) {
            out.perVertexAcceleration.push_back(reg.get<cardillo::C_LinearAcceleration3>(node).value);
        } else {
            out.perVertexAcceleration.emplace_back(Vector3r::Zero());
        }
    }

    out.triangles = surf.triangles;
    return !out.vertices.empty();
}

bool buildPlaneMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& pl = reg.get<cardillo::C_Plane>(e);
    const auto state = cardillo::RigidBody::getState(reg, e);
    const auto& pos = state.position;
    out.center = pos;

    Vector3r n = pl.normal.normalized();
    Vector3r a = pl.up - pl.up.dot(n) * n;
    if (a.norm() < (real_t)1e-12) {
        a = Vector3r::UnitX() - Vector3r::UnitX().dot(n) * n;
    }
    a.normalize();
    const Vector3r b = n.cross(a);

    out.vertices = {pos + (-pl.sizeX) * a + (-pl.sizeY) * b, pos + (pl.sizeX) * a + (-pl.sizeY) * b, pos + (pl.sizeX) * a + (pl.sizeY) * b, pos + (-pl.sizeX) * a + (pl.sizeY) * b};
    out.triangles.emplace_back(0, 1, 2);
    out.triangles.emplace_back(0, 2, 3);

    out.hasUV = true;
    out.uvs = {Eigen::Vector2f(0.f, 0.f), Eigen::Vector2f(static_cast<float>(2 * pl.sizeX), 0.f), Eigen::Vector2f(static_cast<float>(2 * pl.sizeX), static_cast<float>(2 * pl.sizeY)),
               Eigen::Vector2f(0.f, static_cast<float>(2 * pl.sizeY))};
    return true;
}

bool buildCubeMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& c = reg.get<cardillo::C_Cube>(e);
    Quaternion4r qLocal = c.q;
    Matrix33r Rw = out.R * qLocal.toRotationMatrix();
    out.R = Rw;
    out.center = out.center + (Rw * c.center);

    const Vector3r ex = Rw * Vector3r::UnitX();
    const Vector3r ey = Rw * Vector3r::UnitY();
    const Vector3r ez = Rw * Vector3r::UnitZ();
    const Vector3r vx = ex * c.halfExtents.x();
    const Vector3r vy = ey * c.halfExtents.y();
    const Vector3r vz = ez * c.halfExtents.z();
    out.vertices = {out.center - vx - vy - vz, out.center + vx - vy - vz, out.center + vx + vy - vz, out.center - vx + vy - vz,
                    out.center - vx - vy + vz, out.center + vx - vy + vz, out.center + vx + vy + vz, out.center - vx + vy + vz};

    out.triangles = {Eigen::Vector3i(0, 1, 2), Eigen::Vector3i(0, 2, 3), Eigen::Vector3i(4, 6, 5), Eigen::Vector3i(4, 7, 6), Eigen::Vector3i(0, 3, 7), Eigen::Vector3i(0, 7, 4),
                     Eigen::Vector3i(1, 5, 6), Eigen::Vector3i(1, 6, 2), Eigen::Vector3i(0, 4, 5), Eigen::Vector3i(0, 5, 1), Eigen::Vector3i(3, 2, 6), Eigen::Vector3i(3, 6, 7)};

    out.hasUV = true;
    out.uvs.resize(8);
    const float hx = static_cast<float>(c.halfExtents.x());
    const float hy = static_cast<float>(c.halfExtents.y());
    for (int i = 0; i < 8; ++i) {
        float lx = (i == 1 || i == 2 || i == 5 || i == 6) ? hx : -hx;
        float ly = (i == 2 || i == 3 || i == 6 || i == 7) ? hy : -hy;
        out.uvs[(size_t)i] = Eigen::Vector2f(lx + hx, ly + hy);
    }
    return true;
}

bool buildCapsuleMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& c = reg.get<cardillo::C_Capsule>(e);
    std::vector<Vector3r> localV;
    cardillo::io::MeshGenerator::generateCapsuleMesh(8, 3, 1, c.radius, c.halfLength, localV, out.triangles);
    out.vertices.reserve(localV.size());
    for (const auto& v : localV) out.vertices.push_back(out.R * v + out.center);
    return !out.vertices.empty();
}

bool buildCylinderMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& c = reg.get<cardillo::C_Cylinder>(e);
    std::vector<Vector3r> localV;
    cardillo::io::MeshGenerator::generateCylinderMesh(16, c.radius, c.halfLength, localV, out.triangles);
    out.vertices.reserve(localV.size());
    for (const auto& v : localV) out.vertices.push_back(out.R * v + out.center);
    return !out.vertices.empty();
}

bool buildConeMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const auto& c = reg.get<cardillo::C_Cone>(e);
    std::vector<Vector3r> localV;
    cardillo::io::MeshGenerator::generateConeMesh(24, c.radius, c.height, localV, out.triangles);
    out.vertices.reserve(localV.size());
    for (const auto& v : localV) out.vertices.push_back(out.R * v + out.center);
    return !out.vertices.empty();
}

bool buildMeshAssetTriangles(const cardillo::World& sys, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    try {
        const auto& asset = sys.getMeshAsset(e);
        if (!asset.bvh || !asset.bvh->vertices || !asset.bvh->tri_indices) return false;

        const auto& V = *asset.bvh->vertices;
        const auto& F = *asset.bvh->tri_indices;
        out.vertices.reserve(V.size());
        for (const auto& v : V) {
            Vector3r p((real_t)v[0], (real_t)v[1], (real_t)v[2]);
            out.vertices.push_back(out.R * p + out.center);
        }
        out.triangles.reserve(F.size());
        for (const auto& t : F) out.triangles.emplace_back((int)t[0], (int)t[1], (int)t[2]);

        if (asset.hasUV && asset.uvs.size() == out.vertices.size()) {
            out.hasUV = true;
            out.uvs = asset.uvs;
        }
        return true;
    } catch (...) {
        return false;
    }
}

bool buildHeightFieldTriangles(const cardillo::World& sys, entt::entity e, int heightFieldStride, cardillo::io::MeshGenerator::EntityMesh& out) {
    try {
        const auto& asset = sys.getHeightFieldAsset(e);
        if (!asset.hf) return false;

        const auto& hf = *asset.hf;
        const auto& X = hf.getXGrid();
        const auto& Y = hf.getYGrid();
        const auto& H = hf.getHeights();
        const int rows = (int)H.rows();
        const int cols = (int)H.cols();
        if (rows <= 1 || cols <= 1) return false;

        const int s = std::max(1, heightFieldStride);
        std::vector<int> xIdx;
        std::vector<int> yIdx;
        for (int x = 0; x < cols; x += s) xIdx.push_back(x);
        if (!xIdx.empty() && xIdx.back() != cols - 1) xIdx.push_back(cols - 1);
        for (int y = 0; y < rows; y += s) yIdx.push_back(y);
        if (!yIdx.empty() && yIdx.back() != rows - 1) yIdx.push_back(rows - 1);

        const int nxs = (int)xIdx.size();
        const int nys = (int)yIdx.size();
        out.vertices.reserve((size_t)nxs * (size_t)nys);
        out.uvs.reserve((size_t)nxs * (size_t)nys);

        for (int yi = 0; yi < nys; ++yi) {
            int y = yIdx[yi];
            float vf = nys > 1 ? (float)yi / (float)(nys - 1) : 0.f;
            for (int xi = 0; xi < nxs; ++xi) {
                int x = xIdx[xi];
                float uf = nxs > 1 ? (float)xi / (float)(nxs - 1) : 0.f;
                Vector3r pLocal((real_t)X[x], (real_t)Y[y], (real_t)H(y, x));
                out.vertices.push_back(out.R * pLocal + out.center);
                out.uvs.emplace_back(uf, 1.0f - vf);
            }
        }
        out.hasUV = true;

        out.triangles.reserve((size_t)(nxs - 1) * (size_t)(nys - 1) * 2);
        for (int yi = 0; yi < nys - 1; ++yi) {
            for (int xi = 0; xi < nxs - 1; ++xi) {
                int i0 = yi * nxs + xi;
                int i1 = i0 + 1;
                int i2 = i0 + nxs;
                int i3 = i2 + 1;
                out.triangles.emplace_back(i0, i1, i2);
                out.triangles.emplace_back(i1, i3, i2);
            }
        }
        return true;
    } catch (...) {
        return false;
    }
}

bool buildSphereMeshTriangles(const entt::registry& reg, entt::entity e, cardillo::io::MeshGenerator::EntityMesh& out) {
    const real_t radius = reg.get<cardillo::C_Radius>(e).r;
    std::vector<Vector3r> unitVerts;
    std::vector<Eigen::Vector3i> unitTris;
    cardillo::io::MeshGenerator::generateUVSphere(12, 24, unitVerts, unitTris);

    out.vertices.reserve(unitVerts.size());
    out.uvs.reserve(unitVerts.size());
    for (const auto& v : unitVerts) {
        out.vertices.push_back(out.R * (radius * v) + out.center);
        float u = static_cast<float>(std::atan2((double)v.y(), (double)v.x()) / (2.0 * M_PI) + 0.5);
        float t = static_cast<float>(std::acos((double)v.z()) / M_PI);
        out.uvs.emplace_back(u, t);
    }
    out.hasUV = true;
    out.triangles = std::move(unitTris);
    return true;
}

}  // namespace

bool MeshGenerator::buildEntityMesh(const cardillo::World& sys, entt::entity e, int heightFieldStride, EntityMesh& out) {
    const auto& reg = sys.ecs();
    if (!reg.valid(e) || !reg.any_of<cardillo::C_VisualObject>(e)) {
        return false;
    }

    fillCommonEntityData(reg, e, out);

    switch (detectMeshKind(reg, e)) {
        case MeshKind::SoftBody:
            return buildSoftBodyMesh(reg, e, out);
        case MeshKind::Plane:
            return buildPlaneMeshTriangles(reg, e, out);
        case MeshKind::Cube:
            return buildCubeMeshTriangles(reg, e, out);
        case MeshKind::Capsule:
            return buildCapsuleMeshTriangles(reg, e, out);
        case MeshKind::Cylinder:
            return buildCylinderMeshTriangles(reg, e, out);
        case MeshKind::Cone:
            return buildConeMeshTriangles(reg, e, out);
        case MeshKind::MeshAsset:
            return buildMeshAssetTriangles(sys, e, out);
        case MeshKind::HeightField:
            return buildHeightFieldTriangles(sys, e, heightFieldStride, out);
        case MeshKind::Sphere:
            return buildSphereMeshTriangles(reg, e, out);
        case MeshKind::Unsupported:
        default:
            return false;
    }
}

}  // namespace cardillo::io
