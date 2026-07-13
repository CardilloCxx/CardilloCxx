#include "rigid_body_factory.hpp"

#include <coal/shape/geometric_shapes.h>
#include <algorithm>
#include <limits>
#include <optional>
#include <type_traits>

namespace cardillo::physics {
namespace {

bool meshAabbFromAsset(const MeshAsset& asset, Vector3r& center_out, Vector3r& he_out) {
    if (!asset.bvh || !asset.bvh->vertices || asset.bvh->vertices->empty()) return false;

    Vector3r minv = Vector3r::Constant(std::numeric_limits<real_t>::max());
    Vector3r maxv = Vector3r::Constant(std::numeric_limits<real_t>::lowest());
    const Matrix33r Rpa = asset.Rpa;

    for (const auto& v : *asset.bvh->vertices) {
        Vector3r vp((real_t)v[0], (real_t)v[1], (real_t)v[2]);
        vp = Rpa.transpose() * vp;
        minv.x() = std::min(minv.x(), vp.x());
        maxv.x() = std::max(maxv.x(), vp.x());
        minv.y() = std::min(minv.y(), vp.y());
        maxv.y() = std::max(maxv.y(), vp.y());
        minv.z() = std::min(minv.z(), vp.z());
        maxv.z() = std::max(maxv.z(), vp.z());
    }

    center_out = (real_t)0.5 * (maxv + minv);
    he_out = (real_t)0.5 * (maxv - minv);
    return true;
}

inline real_t computeVolume(const RigidShape& shape, const World* system = nullptr) {
    return std::visit(
        [system](auto&& s) -> real_t {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, CubeShape>) {
                return 8.0 * s.halfExtents.x() * s.halfExtents.y() * s.halfExtents.z();
            } else if constexpr (std::is_same_v<T, SphereShape>) {
                return (4.0 / 3.0) * M_PI * s.radius * s.radius * s.radius;
            } else if constexpr (std::is_same_v<T, CylinderShape>) {
                return M_PI * s.radius * s.radius * (2.0 * s.halfLength);
            } else if constexpr (std::is_same_v<T, CapsuleShape>) {
                return M_PI * s.radius * s.radius * (2.0 * s.halfLength + (4.0 / 3.0) * s.radius);
            } else if constexpr (std::is_same_v<T, ConeShape>) {
                return (1.0 / 3.0) * M_PI * s.radius * s.radius * s.height;
            } else if constexpr (std::is_same_v<T, PlaneShape>) {
                return 0.0;
            } else if constexpr (std::is_same_v<T, MeshShape>) {
                if (!system) return 0.0;
                return system->assets().getMesh(s.path, s.scale, true).volume;
            }
        },
        shape);
}

inline Vector3r computeUnitInertia(const RigidShape& shape, const World* system) {
    return std::visit(
        [system](auto&& s) -> Vector3r {
            using T = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<T, CubeShape>) {
                return Vector3r((s.halfExtents.y() * s.halfExtents.y() + s.halfExtents.z() * s.halfExtents.z()) / 3.0,
                                (s.halfExtents.x() * s.halfExtents.x() + s.halfExtents.z() * s.halfExtents.z()) / 3.0,
                                (s.halfExtents.x() * s.halfExtents.x() + s.halfExtents.y() * s.halfExtents.y()) / 3.0);
            } else if constexpr (std::is_same_v<T, SphereShape>) {
                const real_t c = 2.0 / 5.0 * s.radius * s.radius;
                return Vector3r(c, c, c);
            } else if constexpr (std::is_same_v<T, CylinderShape>) {
                const real_t h = 2.0 * s.halfLength;
                const real_t r2 = s.radius * s.radius;
                const real_t Izz = 0.5 * r2;
                const real_t Ixx = (3.0 * r2 + h * h) / 12.0;
                return Vector3r(Ixx, Ixx, Izz);
            } else if constexpr (std::is_same_v<T, CapsuleShape>) {
                coal::Capsule capsule((coal::CoalScalar)s.radius, (coal::CoalScalar)(2.0 * s.halfLength));
                auto Iunit = capsule.computeMomentofInertia();  // convert to unit inertia
                return Vector3r(Iunit(0, 0), Iunit(1, 1), Iunit(2, 2)) / ((real_t)computeVolume(s, system));
            } else if constexpr (std::is_same_v<T, ConeShape>) {
                coal::Cone cone((coal::CoalScalar)s.radius, (coal::CoalScalar)s.height);
                auto Iunit = cone.computeMomentofInertia();  // convert to unit inertia
                return Vector3r(Iunit(0, 0), Iunit(1, 1), Iunit(2, 2)) / ((real_t)computeVolume(s, system));
            } else if constexpr (std::is_same_v<T, PlaneShape>) {
                return Vector3r::Zero();
            } else if constexpr (std::is_same_v<T, MeshShape>) {
                return system->assets().getMesh(s.path, s.scale, true).inertia_diag_unit / ((real_t)computeVolume(s, system));
            }
        },
        shape);
}

inline real_t getMass(const RigidShape& shape, const physics::RigidProps& props, const World* system) {
    // Use explicit mass if provided
    if (props.mass.has_value()) return std::max((real_t)0, props.mass.value());

    // Otherwise, use density if available
    const real_t density = props.density.value_or(0);
    if (density <= 0) return 0;

    const real_t volume = computeVolume(shape, system);
    return density * volume;
}

inline Vector3r getInertia(const RigidShape& shape, real_t mass, const World* system) {
    const Vector3r unitInertia = computeUnitInertia(shape, system);
    return unitInertia * mass;  // scale unit inertia by mass
}

}  // namespace

entt::entity RigidBodyFactory::create(World& system, const physics::RigidShape& shape, const physics::RigidState& state, const physics::RigidProps& props) {
    auto& reg = system.ecs();
    const auto& cfg = system.config();

    const auto e = reg.create();

    // --- Base components ---
    reg.emplace<cardillo::C_Position3>(e, state.position);
    reg.emplace<cardillo::C_Orientation>(e, Quaternion4r(state.orientation).normalized());
    reg.emplace<cardillo::C_LinearVelocity3>(e, state.linearVelocity);
    reg.emplace<cardillo::C_AngularVelocity3>(e, state.angularVelocity);
    reg.emplace<cardillo::C_LinearAcceleration3>(e, Vector3r::Zero());
    reg.emplace<cardillo::C_AngularAcceleration3>(e, Vector3r::Zero());

    const bool wantsColliderVisual = std::holds_alternative<MeshShape>(shape) && std::get<MeshShape>(shape).show_collider;

    if (props.visual || wantsColliderVisual) reg.emplace<cardillo::C_VisualObject>(e);

    if (props.collidable && !cfg.collision_disable_all) reg.emplace<cardillo::C_Collidable>(e);

    const real_t mu = (props.friction >= (real_t)0) ? props.friction : cfg.friction_default_mu;
    reg.emplace<cardillo::C_Friction>(e, mu);

    const real_t restN = (props.restitution_normal >= (real_t)0) ? props.restitution_normal : cfg.restitution_default_normal;
    const real_t restT = (props.restitution_tangential >= (real_t)0) ? props.restitution_tangential : cfg.restitution_default_tangential;
    reg.emplace<cardillo::C_Restitution>(e, std::max((real_t)0, restN), std::max((real_t)0, restT));

    // --- Helper for adding dynamic tags of rigid body ---
    const auto addRigidBody = [&](real_t mass, const Vector3r& inertiaDiag) {
        reg.emplace<cardillo::C_PhysicsObject>(e);
        reg.emplace<cardillo::C_RigidBodyTag>(e);
        reg.emplace<cardillo::C_Mass>(e, mass);
        reg.emplace<cardillo::C_InertiaDiag>(e, inertiaDiag);
    };

    struct Visitor {
        entt::registry& reg;
        World& system;
        entt::entity e;
        const physics::RigidState& state;
        const physics::RigidProps& props;
        const decltype(addRigidBody)& addRigidBodyFn;
        const decltype(cfg)& cfgRef;

        void operator()(const CubeShape& s) const {
            const real_t mass = getMass(s, props, &system);

            if (props.visual) reg.emplace<cardillo::C_CubeVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Cube>(e, Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity());
            reg.emplace<cardillo::C_Cube>(e, Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity());

            if (mass > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
        }

        void operator()(const SphereShape& s) const {
            const real_t mass = getMass(s, props, &system);

            if (props.visual) reg.emplace<cardillo::C_PointVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Sphere>(e);
            reg.emplace<cardillo::C_Radius>(e, s.radius);

            if (mass > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
        }

        void operator()(const CylinderShape& s) const {
            const real_t mass = getMass(s, props, &system);

            if (props.visual) reg.emplace<cardillo::C_CylinderVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Cylinder>(e, s.radius, s.halfLength);
            reg.emplace<cardillo::C_Cylinder>(e, s.radius, s.halfLength);

            if (mass > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
        }

        void operator()(const CapsuleShape& s) const {
            const real_t mass = getMass(s, props, &system);

            if (props.visual) reg.emplace<cardillo::C_CapsuleVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Capsule>(e, s.radius, s.halfLength);
            reg.emplace<cardillo::C_Capsule>(e, s.radius, s.halfLength);

            if (mass > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
        }

        void operator()(const ConeShape& s) const {
            const real_t mass = getMass(s, props, &system);

            if (props.visual) reg.emplace<cardillo::C_ConeVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Cone>(e, s.radius, s.height);
            reg.emplace<cardillo::C_Cone>(e, s.radius, s.height);

            if (mass > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
        }

        void operator()(const PlaneShape& s) const {
            if (props.visual) reg.emplace<cardillo::C_PlaneVisualTag>(e);
            if (props.collidable && !cfgRef.collision_disable_all) reg.emplace<cardillo::C_RB_Plane>(e, s.normal, s.up, s.sizeX, s.sizeY);
            reg.emplace<cardillo::C_Plane>(e, s.normal, s.up, s.sizeX, s.sizeY);
        }

        void operator()(const MeshShape& s) const {
            if (props.visual) reg.emplace<cardillo::C_MeshVisualTag>(e);
            reg.emplace<cardillo::C_Mesh>(e, s.path, s.scale);

            const real_t mass = getMass(s, props, &system);
            const bool dynamic = mass > 0;
            const MeshAsset& asset = system.assets().getMesh(s.path, s.scale, dynamic);

            if (props.collidable && !cfgRef.collision_disable_all) {
                if (s.use_bbox_collider) {
                    Vector3r center = Vector3r::Zero();
                    Vector3r he((real_t)0.05, (real_t)0.05, (real_t)0.05);
                    meshAabbFromAsset(asset, center, he);
                    Quaternion4r q_local(asset.Rpa.transpose());
                    reg.emplace<cardillo::C_RB_Cube>(e, center, he, q_local);

                    if (s.show_collider) {
                        reg.emplace<cardillo::C_Cube>(e, center, he, q_local);
                        reg.emplace<cardillo::C_CubeVisualTag>(e);
                    }
                } else {
                    reg.emplace<cardillo::C_RB_Mesh>(e);
                }
            }

            if (dynamic) {
                const Quaternion4r q_new = state.orientation * Quaternion4r(asset.Rpa);
                const Vector3r pos_new = state.position + (state.orientation * asset.com);
                reg.get<cardillo::C_Position3>(e).value = pos_new;
                reg.get<cardillo::C_Orientation>(e).value = q_new.normalized();

                if (asset.volume > 0) addRigidBodyFn(mass, getInertia(s, mass, &system));
            }
        }
    };

    std::visit(Visitor{reg, system, e, state, props, addRigidBody, cfg}, shape);

    system.markStructureDirty();
    return e;
}

}  // namespace cardillo::physics
