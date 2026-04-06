#pragma once

#include <entt/entt.hpp>
#include <vector>
#include "misc/types.hpp"
#include "../physics/world.hpp"

namespace cardillo::io {

class MeshGenerator {
public:
    struct EntityMesh {
        entt::entity entity{entt::null};
        bool isDynamic{false};
        int entityId{-1};
        real_t beamLengthRatio{(real_t)1};

        std::vector<Vector3r> vertices;
        std::vector<Eigen::Vector3i> triangles;

        bool hasUV{false};
        std::vector<Eigen::Vector2f> uvs;

        bool hasKinematics{false};
        Vector3r center{Vector3r::Zero()};
        Matrix33r R{Matrix33r::Identity()};
        Vector3r vlin{Vector3r::Zero()};
        Vector3r omega{Vector3r::Zero()};
        Vector3r alin{Vector3r::Zero()};
        Vector3r alpha{Vector3r::Zero()};
        real_t entityPressure{(real_t)0};

        std::vector<Vector3r> perVertexVelocity;
        std::vector<Vector3r> perVertexAcceleration;
    };

    // Generate a unit UV-sphere triangulation.
    static void generateUVSphere(int latSegments,
                                 int lonSegments,
                                 std::vector<Vector3r>& outVertices,
                                 std::vector<Eigen::Vector3i>& outTriangles);

    static void generateCapsuleMesh(int segmentsCircumference,
                                    int segmentsHemisphere,
                                    int segmentsCylinder,
                                    real_t radius,
                                    real_t halfLength,
                                    std::vector<Vector3r>& vertices,
                                    std::vector<Eigen::Vector3i>& triangles);

    static void generateCylinderMesh(int segmentsCircumference,
                                     real_t radius,
                                     real_t halfLength,
                                     std::vector<Vector3r>& vertices,
                                     std::vector<Eigen::Vector3i>& triangles);

    static void generateConeMesh(int segmentsCircumference,
                                 real_t radius,
                                 real_t height,
                                 std::vector<Vector3r>& vertices,
                                 std::vector<Eigen::Vector3i>& triangles);

    // Extract a renderable mesh for a visual entity in world coordinates.
    // Returns false if the entity has no supported visual geometry.
    static bool buildEntityMesh(const cardillo::World& sys,
                                entt::entity e,
                                int heightFieldStride,
                                EntityMesh& out);
};

}
