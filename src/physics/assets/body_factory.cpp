#include "body_factory.hpp"

#include <algorithm>
#include <limits>
#include <optional>

#include "../../collision/collision_coal.hpp"

#include <coal/shape/geometric_shapes.h>

#include "../../io/csv_writer.hpp"
#include "../../io/softbody_loader.hpp"
#include "../constraints/constraint_factory.hpp"
#include "../constraints/constraints.hpp"
#include "assets.hpp"

namespace cardillo {
namespace physics {

entt::entity BodyFactory::addPointMass(World& sys, real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius) {
    auto& reg = sys.ecs();
    auto e = reg.create();
    reg.emplace<C_PhysicsObject>(e);
    reg.emplace<C_PointMassTag>(e);
    if (!sys.config().collision_disable_all) reg.emplace<C_Collidable>(e);
    reg.emplace<C_VisualObject>(e);
    reg.emplace<C_PointVisualTag>(e);
    reg.emplace<C_Mass>(e, C_Mass{mass});
    reg.emplace<C_Position3>(e, C_Position3{x0});
    reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{v0});
    reg.emplace<C_LinearAcceleration3>(e, C_LinearAcceleration3{Vector3r::Zero()});
    reg.emplace<C_Radius>(e, C_Radius{radius});
    if (!reg.any_of<C_Friction>(e)) {
        reg.emplace<C_Friction>(e, C_Friction{sys.config().friction_default_mu});
    }
    if (!reg.any_of<C_Restitution>(e)) {
        const auto& cfg = sys.config();
        reg.emplace<C_Restitution>(e, C_Restitution{std::max((real_t)0, cfg.restitution_default_normal), std::max((real_t)0, cfg.restitution_default_tangential)});
    }
    sys.markStructureDirty();
    return e;
}

std::vector<entt::entity> BodyFactory::addSoftBody(World& sys, const std::string& objPath, real_t stiffness, real_t damping, const Vector3r& position, const Quaternion4r& orientation,
                                                   const Vector3r& linearVelocity, const Vector3r& angularVelocity, real_t totalMass, real_t nodeRadius,
                                                   collision::CollisionCoal* collision_mgr) {
    std::vector<entt::entity> nodes;
    io::SoftBodyMesh sb;
    if (!io::load_obj_softbody(objPath, sb)) {
        std::printf("[SoftBody] Failed to load OBJ: %s\n", objPath.c_str());
        return nodes;
    }

    const size_t N = sb.positions.size();
    if (N == 0) return nodes;

    real_t nodeMass = (totalMass > (real_t)0) ? (totalMass / (real_t)N) : (real_t)0.02;

    Matrix33r R = orientation.toRotationMatrix();
    const Vector3r angularVelocityWorld = R * angularVelocity;  // angularVelocity is given in the body-fixed frame
    nodes.reserve(N);
    for (const auto& p0 : sb.positions) {
        Vector3r pw = position + R * p0;
        Vector3r vw = linearVelocity + angularVelocityWorld.cross(pw - position);
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

    if (collision_mgr) {
        for (const auto& n1 : nodes) {
            for (const auto& n2 : nodes) {
                if (n1 < n2) {
                    collision_mgr->disablePair(n1, n2);
                }
            }
        }
    }

    if (!sb.triangles.empty()) {
        entt::entity surf = sys.ecs().create();
        sys.ecs().emplace<C_VisualObject>(surf);
        sys.ecs().emplace<C_SoftBodyVisualTag>(surf);
        C_SoftBodySurface surfComp;
        surfComp.triangles = sb.triangles;
        surfComp.nodes = nodes;
        sys.ecs().emplace<C_SoftBodySurface>(surf, std::move(surfComp));
    }

    sys.markStructureDirty();
    return nodes;
}

}  // namespace physics
}  // namespace cardillo
