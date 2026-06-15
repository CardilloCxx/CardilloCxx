#include "body_factory.hpp"

#include <algorithm>
#include <limits>
#include <optional>

#include "../../collision/collision_coal.hpp"

#include <coal/shape/geometric_shapes.h>

#include "../../io/csv_writer.hpp"
#include "../../io/heightmap_loader.hpp"
#include "../../io/softbody_loader.hpp"
#include "../constraints/constraint_factory.hpp"
#include "../constraints/constraints.hpp"
#include "assets.hpp"

namespace cardillo {
namespace physics {

entt::entity BodyFactory::addPointMass(World& sys, real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius) {
    auto& reg = sys.ecs();
    auto e = reg.create();
    reg.emplace<cardillo::C_PhysicsObject>(e);
    reg.emplace<cardillo::C_PointMassTag>(e);
    if (!sys.config().collision_disable_all) reg.emplace<cardillo::C_Collidable>(e);
    reg.emplace<cardillo::C_VisualObject>(e);
    reg.emplace<cardillo::C_PointVisualTag>(e);
    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
    reg.emplace<cardillo::C_Position3>(e, cardillo::C_Position3{x0});
    reg.emplace<cardillo::C_LinearVelocity3>(e, cardillo::C_LinearVelocity3{v0});
    reg.emplace<cardillo::C_LinearAcceleration3>(e, cardillo::C_LinearAcceleration3{Vector3r::Zero()});
    reg.emplace<cardillo::C_Radius>(e, cardillo::C_Radius{radius});
    if (!reg.any_of<cardillo::C_Friction>(e)) {
        reg.emplace<cardillo::C_Friction>(e, cardillo::C_Friction{sys.config().friction_default_mu});
    }
    sys.markStructureDirty();
    return e;
}

entt::entity BodyFactory::addObstacleHeightField(World& sys, const Vector3r& position, const Quaternion4r& orientation, const std::string& exrPath, real_t x_dim, real_t y_dim, real_t z_scale,
                                            real_t min_height) {
    auto& reg = sys.ecs();
    auto e = reg.create();
    if (!sys.config().collision_disable_all) reg.emplace<cardillo::C_Collidable>(e);
    reg.emplace<cardillo::C_VisualObject>(e);
    reg.emplace<cardillo::C_Position3>(e, cardillo::C_Position3{position});
    reg.emplace<cardillo::C_Orientation>(e, cardillo::C_Orientation{orientation});
    reg.emplace<cardillo::C_HeightFieldVisualTag>(e);
    reg.emplace<cardillo::C_HeightField>(e, cardillo::C_HeightField{exrPath, x_dim, y_dim, z_scale, min_height});
    if (!sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_HeightField>(e);
    if (!reg.any_of<cardillo::C_Friction>(e)) {
        reg.emplace<cardillo::C_Friction>(e, cardillo::C_Friction{sys.config().friction_default_mu});
    }
    sys.markStructureDirty();
    (void)sys.getHeightFieldAsset(e);
    return e;
}

std::vector<entt::entity> BodyFactory::addSoftBody(World& sys, const std::string& objPath, real_t stiffness, real_t damping, const Vector3r& position, const Quaternion4r& orientation,
                                                   const Vector3r& linearVelocity, const Vector3r& angularVelocity, real_t totalMass, real_t nodeRadius,
                                                   cardillo::collision::CollisionCoal* collision_mgr) {
    std::vector<entt::entity> nodes;
    cardillo::io::SoftBodyMesh sb;
    if (!cardillo::io::load_obj_softbody(objPath, sb)) {
        std::printf("[SoftBody] Failed to load OBJ: %s\n", objPath.c_str());
        return nodes;
    }

    const size_t N = sb.positions.size();
    if (N == 0) return nodes;

    real_t nodeMass = (totalMass > (real_t)0) ? (totalMass / (real_t)N) : (real_t)0.02;

    Matrix33r R = orientation.toRotationMatrix();
    nodes.reserve(N);
    for (const auto& p0 : sb.positions) {
        Vector3r pw = position + R * p0;
        Vector3r vw = linearVelocity + angularVelocity.cross(pw - position);
        entt::entity id = addPointMass(sys, nodeMass, pw, vw, nodeRadius);
        nodes.push_back(id);
    }

    for (const auto& e : sb.edges) {
        int i = e.first;
        int j = e.second;
        if (i >= 0 && j >= 0 && (size_t)i < nodes.size() && (size_t)j < nodes.size()) {
            entt::entity A = nodes[(size_t)i];
            entt::entity B = nodes[(size_t)j];
            ConstraintFactory::addLinearDistanceConstraint(sys, A, B, Vector3r::Zero(), Vector3r::Zero(), stiffness, damping);
        }
    }

    for (const auto& n1 : nodes) {
        for (const auto& n2 : nodes) {
            if (n1 < n2) {
                collision_mgr->disablePair(n1, n2);
            }
        }
    }

    if (!sb.triangles.empty()) {
        entt::entity surf = sys.ecs().create();
        sys.ecs().emplace<cardillo::C_VisualObject>(surf);
        sys.ecs().emplace<cardillo::C_SoftBodyVisualTag>(surf);
        cardillo::C_SoftBodySurface surfComp;
        surfComp.triangles = sb.triangles;
        surfComp.nodes = nodes;
        sys.ecs().emplace<cardillo::C_SoftBodySurface>(surf, std::move(surfComp));
    }

    sys.markStructureDirty();
    return nodes;
}

}  // namespace physics
}  // namespace cardillo
