#include "body_factory.hpp"

#include <algorithm>
#include <limits>
#include <optional>

#include "../../collision/collision_coal.hpp"

#include <coal/shape/geometric_shapes.h>

#include "assets.hpp"
#include "../constraints/constraints.hpp"
#include "../constraints/constraint_factory.hpp"
#include "../../io/csv_writer.hpp"
#include "../../io/heightmap_loader.hpp"
#include "../../io/softbody_loader.hpp"

namespace cardillo {
namespace physics {

namespace {

bool meshAabbFromAsset(const MeshAsset& asset, Vector3r& center_out, Vector3r& he_out) {
    if (!asset.bvh || !asset.bvh->vertices || asset.bvh->vertices->empty()) return false;
    Vector3r minv = Vector3r::Constant(std::numeric_limits<real_t>::max());
    Vector3r maxv = Vector3r::Constant(std::numeric_limits<real_t>::lowest());
    const Matrix33r Rpa = asset.Rpa;
    for (const auto& v : *asset.bvh->vertices) {
        Vector3r vp((real_t)v[0], (real_t)v[1], (real_t)v[2]);
        vp = Rpa.transpose() * vp;
        minv.x() = std::min(minv.x(), vp.x()); maxv.x() = std::max(maxv.x(), vp.x());
        minv.y() = std::min(minv.y(), vp.y()); maxv.y() = std::max(maxv.y(), vp.y());
        minv.z() = std::min(minv.z(), vp.z()); maxv.z() = std::max(maxv.z(), vp.z());
    }
    center_out = (real_t)0.5 * (maxv + minv);
    he_out = (real_t)0.5 * (maxv - minv);
    return true;
}

Vector3r boxInertiaDiag(real_t m, const Vector3r& he) {
    Vector3r I;
    I.x() = (real_t)1.0 / 3.0 * m * (he.y() * he.y() + he.z() * he.z());
    I.y() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.z() * he.z());
    I.z() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.y() * he.y());
    return I;
}

Vector3r sphereInertiaDiag(real_t m, real_t r) {
    const real_t c = (real_t)2.0 / 5.0 * m * r * r;
    return Vector3r(c, c, c);
}

Vector3r capsuleInertiaDiag(real_t m, real_t radius, real_t halfLength) {
    coal::Capsule capsule((coal::CoalScalar)radius, (coal::CoalScalar)(halfLength * 2));
    const real_t volume = static_cast<real_t>(capsule.computeVolume());
    if (volume <= (real_t)0) return Vector3r::Zero();
    const auto Iunit = capsule.computeMomentofInertia();
    const real_t scale = m / volume;
    return Vector3r(static_cast<real_t>(Iunit(0, 0)),
                    static_cast<real_t>(Iunit(1, 1)),
                    static_cast<real_t>(Iunit(2, 2))) * scale;
}

Vector3r cylinderInertiaDiag(real_t m, real_t radius, real_t halfLength) {
    const real_t h = halfLength * (real_t)2.0;
    const real_t r2 = radius * radius;
    const real_t Izz = (real_t)0.5 * m * r2;
    const real_t Ixx = (real_t)1.0 / 12.0 * m * ((real_t)3.0 * r2 + h * h);
    return Vector3r(Ixx, Ixx, Izz);
}

Vector3r coneInertiaDiag(real_t m, real_t radius, real_t height) {
    coal::Cone cone((coal::CoalScalar)radius, (coal::CoalScalar)height);
    const real_t volume = static_cast<real_t>(cone.computeVolume());
    if (volume <= (real_t)0) return Vector3r::Zero();
    const auto Iunit = cone.computeMomentofInertia();
    const real_t scale = m / volume;
    return Vector3r(static_cast<real_t>(Iunit(0, 0)),
                    static_cast<real_t>(Iunit(1, 1)),
                    static_cast<real_t>(Iunit(2, 2))) * scale;
}

struct BeamSample {
    Vector3r position;
    Vector3r tangent;
    Vector3r normal;
    Vector3r binormal;
    real_t segLen;
};

Matrix33r makeFrameFromTangentLocal(const Vector3r& tangent) {
    Vector3r T = tangent.normalized();
    Vector3r up = std::abs(T.z()) < (real_t)0.9 ? Vector3r(0, 0, 1) : Vector3r(0, 1, 0);
    Vector3r N = up.cross(T).normalized();
    if (N.squaredNorm() < (real_t)1e-12) {
        up = Vector3r(1, 0, 0);
        N = up.cross(T).normalized();
    }
    Vector3r B = T.cross(N).normalized();
    Matrix33r M;
    M.col(0) = T;
    M.col(1) = N;
    M.col(2) = B;
    return M;
}

std::pair<entt::entity, entt::entity> buildBeamFromSamples(World& sys,
                                                            const std::vector<BeamSample>& samples,
                                                            bool loop,
                                                            const World::BeamCrossSection& section,
                                                            const World::BeamSpringParams& springs,
                                                            const World::RigidState& stateDefaults,
                                                            const World::RigidProps& propsDefaults,
                                                            const Vector3r& splineCOMWorld,
                                                            cardillo::collision::CollisionCoal* collision_mgr) {
    if (samples.empty()) return {entt::null, entt::null};

    real_t totalLen = (real_t)0;
    for (const auto& s : samples) totalLen += s.segLen;
    if (totalLen <= (real_t)0) totalLen = (real_t)1;

    Matrix33r Rshape = Matrix33r::Identity();
    if (section.type == World::BeamBodyType::Capsule || section.type == World::BeamBodyType::Cylinder) {
        Rshape = Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitX()).toRotationMatrix();
    }

    entt::entity root = entt::null;
    entt::entity prev = entt::null;
    entt::entity end = entt::null;

    Matrix33r Rbody = stateDefaults.orientation.toRotationMatrix();
    Vector3r worldCOM = splineCOMWorld + stateDefaults.position;
    Vector3r v_body_world = Rbody * stateDefaults.linearVelocity;
    Vector3r w_body_world = Rbody * stateDefaults.angularVelocity;

    Quaternion4r q_prev = Quaternion4r::Identity();

    for (const auto& s : samples) {
        const real_t segLen = s.segLen;
        World::RigidShape shape;
        if (section.type == World::BeamBodyType::Cube) {
            shape = CubeShape(Vector3r(segLen * (real_t)0.5, section.width * (real_t)0.5, section.height * (real_t)0.5));
        } else if (section.type == World::BeamBodyType::Cylinder) {
            real_t r = std::min(section.width, section.height) * (real_t)0.5;
            shape = CylinderShape(r, segLen * (real_t)0.5);
        } else {
            real_t r = std::min(section.width, section.height) * (real_t)0.5;
            shape = CapsuleShape(r, segLen * (real_t)0.5);
        }

        World::RigidProps segProps = propsDefaults;
        real_t massPerSegment = 0;
        if (propsDefaults.mass.has_value()) massPerSegment = *propsDefaults.mass * (segLen / totalLen);
        else if (propsDefaults.density.has_value()) massPerSegment = *propsDefaults.density * (section.area() * segLen);
        segProps.mass = (massPerSegment > 0 ? std::optional<real_t>(massPerSegment) : std::nullopt);

        Matrix33r Rlocal;
        if (s.normal.squaredNorm() > (real_t)0 && s.binormal.squaredNorm() > (real_t)0) {
            Rlocal.col(0) = s.tangent.normalized();
            Rlocal.col(1) = s.normal.normalized();
            Rlocal.col(2) = s.binormal.normalized();
        } else {
            Rlocal = makeFrameFromTangentLocal(s.tangent);
        }

        Matrix33r Rworld = Rbody * Rlocal * Rshape;
        Quaternion4r qworld(Rworld);
        qworld.normalize();
        qworld = World::alignQuaternionTo(qworld, q_prev);
        q_prev = qworld;

        Vector3r worldPos = splineCOMWorld + stateDefaults.position + Rbody * (s.position - splineCOMWorld);
        Vector3r v_world = v_body_world + w_body_world.cross(worldPos - worldCOM);

        World::RigidState segState;
        segState.position = worldPos;
        segState.orientation = qworld;
        segState.linearVelocity = v_world;
        segState.angularVelocity = Rlocal.transpose() * stateDefaults.angularVelocity;

        entt::entity cur = BodyFactory::addRigidBody(sys, shape, segState, segProps);

        cardillo::C_BeamElement be_cur;
        be_cur.l0 = segLen;
        be_cur.l = segLen;
        if (prev != entt::null) {
            be_cur.prev = prev;
            if (!sys.ecs().any_of<cardillo::C_BeamElement>(prev)) {
                cardillo::C_BeamElement be_prev;
                be_prev.l0 = segLen;
                be_prev.l = segLen;
                be_prev.next = cur;
                sys.ecs().emplace<cardillo::C_BeamElement>(prev, be_prev);
            } else {
                sys.ecs().get<cardillo::C_BeamElement>(prev).next = cur;
            }
        }
        sys.ecs().emplace<cardillo::C_BeamElement>(cur, be_cur);

        if (prev != entt::null) {
            ConstraintFactory::addBeamConstraint(sys, prev, cur, springs, section);
            if (collision_mgr) collision_mgr->disablePair(prev, cur);
        }
        if (root == entt::null) root = cur;
        prev = cur;
        end = cur;
    }

        if (loop && root != entt::null && end != entt::null && end != root) {
        ConstraintFactory::addBeamConstraint(sys, end, root, springs, section);
        if (collision_mgr) collision_mgr->disablePair(end, root);
        if (sys.ecs().any_of<cardillo::C_BeamElement>(end)) sys.ecs().get<cardillo::C_BeamElement>(end).next = root;
        if (sys.ecs().any_of<cardillo::C_BeamElement>(root)) sys.ecs().get<cardillo::C_BeamElement>(root).prev = end;
    }

    return {root, end};
}

} // namespace

entt::entity BodyFactory::addRigidBody(World& sys,
                                       const World::RigidShape& shape,
                                       const World::RigidState& state,
                                       const World::RigidProps& props) {
    auto& reg = sys.ecs();
    auto e = reg.create();

    reg.emplace<cardillo::C_Position3>(e, cardillo::C_Position3{state.position});
    reg.emplace<cardillo::C_Orientation>(e, cardillo::C_Orientation{Quaternion4r(state.orientation).normalized()});
    reg.emplace<cardillo::C_LinearVelocity3>(e, cardillo::C_LinearVelocity3{state.linearVelocity});
    reg.emplace<cardillo::C_AngularVelocity3>(e, cardillo::C_AngularVelocity3{state.angularVelocity});

    bool wantsColliderVisual = std::holds_alternative<MeshShape>(shape)
        && std::get<MeshShape>(shape).show_collider;
    if (props.visual || wantsColliderVisual) reg.emplace<cardillo::C_VisualObject>(e);
    if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_Collidable>(e);

    std::optional<real_t> massOpt = props.mass;
    real_t computedMass = (real_t)0;
    auto computeVolumeCube = [](const Vector3r& he) { return (real_t)8 * he.x() * he.y() * he.z(); };
    auto computeVolumeCapsule = [](real_t r, real_t h) { return (real_t)M_PI * r * r * (2 * h + (real_t)4.0 / 3.0 * r); };
    auto computeVolumeCylinder = [](real_t r, real_t h) { return (real_t)M_PI * r * r * (2 * h); };
    auto computeVolumeSphere = [](real_t r) { return (real_t)4.0 / 3.0 * (real_t)M_PI * r * r * r; };
    auto computeVolumeCone = [](real_t r, real_t h) { return (real_t)1.0 / 3.0 * (real_t)M_PI * r * r * h; };
    real_t densityUsed = props.density.value_or((real_t)0);

    if (!massOpt.has_value() && props.density.has_value()) {
        std::visit([&](auto&& s) {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, CubeShape>) computedMass = densityUsed * computeVolumeCube(s.halfExtents);
            else if constexpr (std::is_same_v<T, CapsuleShape>) computedMass = densityUsed * computeVolumeCapsule(s.radius, s.halfLength);
            else if constexpr (std::is_same_v<T, CylinderShape>) computedMass = densityUsed * computeVolumeCylinder(s.radius, s.halfLength);
            else if constexpr (std::is_same_v<T, SphereShape>) computedMass = densityUsed * computeVolumeSphere(s.radius);
            else if constexpr (std::is_same_v<T, ConeShape>) computedMass = densityUsed * computeVolumeCone(s.radius, s.height);
            else if constexpr (std::is_same_v<T, PlaneShape>) computedMass = 0;
            else if constexpr (std::is_same_v<T, MeshShape>) {
                const ::cardillo::MeshAsset& ma = sys.assets().getMesh(s.path, s.scale, true);
                if (ma.volume > (real_t)0) computedMass = densityUsed * ma.volume;
            }
        }, shape);
        if (computedMass > (real_t)0) massOpt = computedMass;
    }

    const real_t mass = std::max((real_t)0, massOpt.value_or((real_t)0));

    real_t mu = (props.friction >= (real_t)0) ? props.friction : sys.config().friction_default_mu;
    reg.emplace<cardillo::C_Friction>(e, cardillo::C_Friction{mu});

        std::visit([&](auto&& s) {
            using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, CubeShape>) {
            if (props.visual) reg.emplace<cardillo::C_CubeVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Cube>(e, cardillo::C_RB_Cube{Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity()});
            reg.emplace<cardillo::C_Cube>(e, cardillo::C_Cube{Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity()});
            if (mass > 0) {
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{boxInertiaDiag(mass, s.halfExtents)});
            }
        } else if constexpr (std::is_same_v<T, SphereShape>) {
            if (props.visual) reg.emplace<cardillo::C_PointVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Sphere>(e);
            reg.emplace<cardillo::C_Radius>(e, cardillo::C_Radius{s.radius});
            if (mass > 0) {
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{sphereInertiaDiag(mass, s.radius)});
            }
        } else if constexpr (std::is_same_v<T, ConeShape>) {
            if (props.visual) reg.emplace<cardillo::C_ConeVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Cone>(e, cardillo::C_RB_Cone{s.radius, s.height});
            reg.emplace<cardillo::C_Cone>(e, cardillo::C_Cone{s.radius, s.height});
            if (mass > 0) {
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{coneInertiaDiag(mass, s.radius, s.height)});
            }
        } else if constexpr (std::is_same_v<T, CapsuleShape>) {
            if (props.visual) reg.emplace<cardillo::C_CapsuleVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Capsule>(e, cardillo::C_RB_Capsule{s.radius, s.halfLength});
            reg.emplace<cardillo::C_Capsule>(e, cardillo::C_Capsule{s.radius, s.halfLength});
            if (mass > 0) {
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{capsuleInertiaDiag(mass, s.radius, s.halfLength)});
            }
        } else if constexpr (std::is_same_v<T, CylinderShape>) {
            if (props.visual) reg.emplace<cardillo::C_CylinderVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Cylinder>(e, cardillo::C_RB_Cylinder{s.radius, s.halfLength});
            reg.emplace<cardillo::C_Cylinder>(e, cardillo::C_Cylinder{s.radius, s.halfLength});
            if (mass > 0) {
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{cylinderInertiaDiag(mass, s.radius, s.halfLength)});
            }
        } else if constexpr (std::is_same_v<T, PlaneShape>) {
            if (props.visual) reg.emplace<cardillo::C_PlaneVisualTag>(e);
            if (props.collidable && !sys.config().collision_disable_all) reg.emplace<cardillo::C_RB_Plane>(e, cardillo::C_RB_Plane{s.normal, s.up, s.sizeX, s.sizeY});
            reg.emplace<cardillo::C_Plane>(e, cardillo::C_Plane{s.normal, s.up, s.sizeX, s.sizeY});
        } else if constexpr (std::is_same_v<T, MeshShape>) {
            if (props.visual) reg.emplace<cardillo::C_MeshVisualTag>(e);
            reg.emplace<cardillo::C_Mesh>(e, cardillo::C_Mesh{s.path, s.scale});

            const bool dynamic = mass > 0;
            const ::cardillo::MeshAsset& asset = sys.assets().getMesh(s.path, s.scale, dynamic);
            if (props.collidable && !sys.config().collision_disable_all) {
                if (s.use_bbox_collider) {
                    Vector3r center = Vector3r::Zero();
                    Vector3r he = Vector3r((real_t)0.05, (real_t)0.05, (real_t)0.05);
                    meshAabbFromAsset(asset, center, he);
                    Quaternion4r q_local(asset.Rpa.transpose());
                    reg.emplace<cardillo::C_RB_Cube>(e, cardillo::C_RB_Cube{center, he, q_local});
                    if (s.show_collider) {
                        reg.emplace<cardillo::C_Cube>(e, cardillo::C_Cube{center, he, q_local});
                        reg.emplace<cardillo::C_CubeVisualTag>(e);
                    }
                } else {
                    reg.emplace<cardillo::C_RB_Mesh>(e);
                }
            }

            if (dynamic) {
                Quaternion4r q_rpa(asset.Rpa);
                Quaternion4r q_new = state.orientation * q_rpa;
                Vector3r pos_new = state.position + (state.orientation * asset.com);
                reg.get<cardillo::C_Position3>(e).value = pos_new;
                reg.get<cardillo::C_Orientation>(e).value = Quaternion4r(q_new).normalized();
                reg.emplace<cardillo::C_PhysicsObject>(e);
                reg.emplace<cardillo::C_RigidBodyTag>(e);
                reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                if (asset.volume > (real_t)0) {
                    const real_t rho = mass / asset.volume;
                    Vector3r Idiag = rho * asset.inertia_diag_unit.cwiseMax(Vector3r::Zero());
                    reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{Idiag});
                }
            }
        }
    }, shape);

    sys.markStructureDirty();
    return e;
}

entt::entity BodyFactory::addStaticBody(World& sys,
                                        const World::RigidShape& shape,
                                        const World::RigidState& state) {
    World::RigidProps staticProps;
    return addRigidBody(sys, shape, state, staticProps);
}

index_t BodyFactory::addPointMass(World& sys,
                                  real_t mass,
                                  const Vector3r& x0,
                                  const Vector3r& v0,
                                  real_t radius) {
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
    reg.emplace<cardillo::C_Radius>(e, cardillo::C_Radius{radius});
    if (!reg.any_of<cardillo::C_Friction>(e)) {
        reg.emplace<cardillo::C_Friction>(e, cardillo::C_Friction{sys.config().friction_default_mu});
    }
    sys.markStructureDirty();
    return static_cast<index_t>(entt::to_integral(e));
}

index_t BodyFactory::addObstacleHeightField(World& sys,
                                            const Vector3r& position,
                                            const Quaternion4r& orientation,
                                            const std::string& exrPath,
                                            real_t x_dim,
                                            real_t y_dim,
                                            real_t z_scale,
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
    return static_cast<index_t>(entt::to_integral(e));
}

std::vector<entt::entity> BodyFactory::addSoftBody(World& sys,
                                                    const std::string& objPath,
                                                    real_t stiffness,
                                                    real_t damping,
                                                    const Vector3r& position,
                                                    const Quaternion4r& orientation,
                                                    const Vector3r& linearVelocity,
                                                    const Vector3r& angularVelocity,
                                                    real_t totalMass) {
    std::vector<entt::entity> nodes;
    cardillo::io::SoftBodyMesh sb;
    if (!cardillo::io::load_obj_softbody(objPath, sb)) {
        std::printf("[SoftBody] Failed to load OBJ: %s\n", objPath.c_str());
        return nodes;
    }

    const size_t N = sb.positions.size();
    if (N == 0) return nodes;

    real_t nodeMass = (totalMass > (real_t)0) ? (totalMass / (real_t)N) : (real_t)0.02;
    const real_t nodeRadius = (real_t)0.02;

    Matrix33r R = orientation.toRotationMatrix();
    nodes.reserve(N);
    for (const auto& p0 : sb.positions) {
        Vector3r pw = position + R * p0;
        Vector3r vw = linearVelocity + angularVelocity.cross(pw - position);
        index_t id = addPointMass(sys, nodeMass, pw, vw, nodeRadius);
        nodes.push_back(entt::entity(static_cast<uint32_t>(id)));
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

std::pair<entt::entity, entt::entity> BodyFactory::createBeam(World& sys,
                                                               const misc::SplinePattern& spline,
                                                               const World::BeamCrossSection& section,
                                                               const World::BeamSpringParams& springs,
                                                               const World::RigidState& stateDefaults,
                                                               const World::RigidProps& propsDefaults,
                                                               size_t segments,
                                                               cardillo::collision::CollisionCoal* collision_mgr) {
    const real_t totalLen = spline.totalLength();
    const real_t segLen = totalLen / (real_t)segments;
    const bool endsOnSpline = true;

    std::vector<BeamSample> samples;
    samples.reserve(segments);

    if (endsOnSpline) {
        const real_t minSegLen = (real_t)1e-8;
        const size_t segCount = segments;
        for (size_t i = 0; i < segCount; ++i) {
            real_t alpha0 = (real_t)i / (real_t)segments;
            real_t alpha1 = (real_t)(i + 1) / (real_t)segments;
            if (spline.isLoop() && i + 1 == segCount) alpha1 = (real_t)0;

            misc::SplineSample si0 = spline.sample(alpha0);
            misc::SplineSample si1 = spline.sample(alpha1);
            Vector3r midPos = (si0.position + si1.position) * (real_t)0.5;
            real_t local_segLen = (si1.position - si0.position).norm();
            if (local_segLen <= minSegLen) continue;
            Vector3r midTangent = (si1.position - si0.position) / local_segLen;
            Vector3r midNormal = (si0.normal + si1.normal).normalized();
            Vector3r midBinormal = (si0.binormal + si1.binormal).normalized();
            samples.push_back(BeamSample{midPos, midTangent, midNormal, midBinormal, local_segLen});
        }
    } else {
        for (size_t i = 0; i <= segments; ++i) {
            real_t alpha = (real_t)i / (real_t)segments;
            misc::SplineSample si = spline.sample(alpha);
            samples.push_back(BeamSample{si.position, si.tangent, si.normal, si.binormal, segLen});
        }
    }

    Vector3r splineCOMWorld = spline.centerOfMass();
    return buildBeamFromSamples(sys, samples, spline.isLoop(), section, springs, stateDefaults, propsDefaults, splineCOMWorld, collision_mgr);
}

std::pair<entt::entity, entt::entity> BodyFactory::createBeams(World& sys,
                                                               const std::vector<const misc::SplinePattern*>& splines,
                                                               const World::BeamCrossSection& section,
                                                               const World::BeamSpringParams& springs,
                                                               const World::RigidState& stateDefaults,
                                                               const World::RigidProps& propsDefaults,
                                                               size_t segmentsPerSpline,
                                                               cardillo::collision::CollisionCoal* collision_mgr) {
    real_t totalLen = 0;
    for (const auto* sp : splines) {
        if (sp) totalLen += sp->totalLength();
    }

    entt::entity first = entt::null;
    entt::entity second = entt::null;
    entt::entity prevEnd = entt::null;
    for (size_t i = 0; i < splines.size(); ++i) {
        auto pair = createBeam(sys,
                       *splines[i],
                       section,
                       springs,
                       stateDefaults,
                       propsDefaults,
                       (size_t)(segmentsPerSpline * (splines[i]->totalLength() / totalLen)),
                       collision_mgr);
        if (first == entt::null) first = pair.first;
        if (prevEnd != entt::null && pair.first != entt::null) {
            if (sys.ecs().any_of<cardillo::C_Orientation>(prevEnd) && sys.ecs().any_of<cardillo::C_Orientation>(pair.first)) {
                auto& qNext = sys.ecs().get<cardillo::C_Orientation>(pair.first).value;
                const auto& qPrev = sys.ecs().get<cardillo::C_Orientation>(prevEnd).value;
                qNext = World::alignQuaternionTo(qNext, qPrev);
            }
            ConstraintFactory::addRigidConstraint(sys, prevEnd, pair.first);
            if (collision_mgr) collision_mgr->disablePair(prevEnd, pair.first);
        }
        prevEnd = pair.second;
    }
    second = prevEnd;
    return {first, second};
}

} // namespace physics
} // namespace cardillo
