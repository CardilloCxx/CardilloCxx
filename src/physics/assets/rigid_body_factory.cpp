#include "rigid_body_factory.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <type_traits>
#include <coal/shape/geometric_shapes.h>

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

inline Vector3r boxInertiaDiag(real_t m, const Vector3r& he) {
    Vector3r I;
    I.x() = (real_t)1.0 / 3.0 * m * (he.y() * he.y() + he.z() * he.z());
    I.y() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.z() * he.z());
    I.z() = (real_t)1.0 / 3.0 * m * (he.x() * he.x() + he.y() * he.y());
    return I;
}

inline Vector3r sphereInertiaDiag(real_t m, real_t r) {
    const real_t c = (real_t)2.0 / 5.0 * m * r * r;
    return Vector3r(c, c, c);
}

inline Vector3r capsuleInertiaDiag(real_t m, real_t radius, real_t halfLength) {
    coal::Capsule capsule((coal::CoalScalar)radius, (coal::CoalScalar)(halfLength * 2));
    const real_t volume = static_cast<real_t>(capsule.computeVolume());
    if (volume <= (real_t)0) return Vector3r::Zero();

    const auto Iunit = capsule.computeMomentofInertia();
    const real_t scale = m / volume;
    return Vector3r(static_cast<real_t>(Iunit(0, 0)),
                    static_cast<real_t>(Iunit(1, 1)),
                    static_cast<real_t>(Iunit(2, 2))) *
           scale;
}

inline Vector3r cylinderInertiaDiag(real_t m, real_t radius, real_t halfLength) {
    const real_t h = halfLength * (real_t)2.0;
    const real_t r2 = radius * radius;
    const real_t Izz = (real_t)0.5 * m * r2;
    const real_t Ixx = (real_t)1.0 / 12.0 * m * ((real_t)3.0 * r2 + h * h);
    const real_t Iyy = Ixx;
    return Vector3r(Ixx, Iyy, Izz);
}

inline Vector3r coneInertiaDiag(real_t m, real_t radius, real_t height) {
    coal::Cone cone((coal::CoalScalar)radius, (coal::CoalScalar)height);
    const real_t volume = static_cast<real_t>(cone.computeVolume());
    if (volume <= (real_t)0) return Vector3r::Zero();

    const auto Iunit = cone.computeMomentofInertia();
    const real_t scale = m / volume;
    return Vector3r(static_cast<real_t>(Iunit(0, 0)),
                    static_cast<real_t>(Iunit(1, 1)),
                    static_cast<real_t>(Iunit(2, 2))) *
           scale;
}

} // namespace

entt::entity RigidBodyFactory::create(World& system,
                                      const World::RigidShape& shape,
                                      const World::RigidState& state,
                                      const World::RigidProps& props) {
    auto& reg = system.ecs();
    const auto& cfg = system.config();

    const auto e = reg.create();

    reg.emplace<cardillo::C_Position3>(e, cardillo::C_Position3{state.position});
    reg.emplace<cardillo::C_Orientation>(e, cardillo::C_Orientation{Quaternion4r(state.orientation).normalized()});
    reg.emplace<cardillo::C_LinearVelocity3>(e, cardillo::C_LinearVelocity3{state.linearVelocity});
    reg.emplace<cardillo::C_AngularVelocity3>(e, cardillo::C_AngularVelocity3{state.angularVelocity});

    const bool wantsColliderVisual = std::holds_alternative<MeshShape>(shape) &&
                                     std::get<MeshShape>(shape).show_collider;
    if (props.visual || wantsColliderVisual) reg.emplace<cardillo::C_VisualObject>(e);
    if (props.collidable && !cfg.collision_disable_all) reg.emplace<cardillo::C_Collidable>(e);

    std::optional<real_t> massOpt = props.mass;
    real_t computedMass = (real_t)0;
    const real_t densityUsed = props.density.value_or((real_t)0);

    const auto computeVolumeCube = [](const Vector3r& he) { return (real_t)8 * he.x() * he.y() * he.z(); };
    const auto computeVolumeCapsule = [](real_t r, real_t h) {
        return (real_t)M_PI * r * r * ((real_t)2.0 * h + (real_t)4.0 / (real_t)3.0 * r);
    };
    const auto computeVolumeCylinder = [](real_t r, real_t h) { return (real_t)M_PI * r * r * ((real_t)2.0 * h); };
    const auto computeVolumeSphere = [](real_t r) { return (real_t)4.0 / (real_t)3.0 * (real_t)M_PI * r * r * r; };
    const auto computeVolumeCone = [](real_t r, real_t h) { return (real_t)1.0 / (real_t)3.0 * (real_t)M_PI * r * r * h; };

    if (!massOpt.has_value() && props.density.has_value()) {
        std::visit(
            [&](auto&& s) {
                using T = std::decay_t<decltype(s)>;
                if constexpr (std::is_same_v<T, CubeShape>) {
                    computedMass = densityUsed * computeVolumeCube(s.halfExtents);
                } else if constexpr (std::is_same_v<T, CapsuleShape>) {
                    computedMass = densityUsed * computeVolumeCapsule(s.radius, s.halfLength);
                } else if constexpr (std::is_same_v<T, CylinderShape>) {
                    computedMass = densityUsed * computeVolumeCylinder(s.radius, s.halfLength);
                } else if constexpr (std::is_same_v<T, SphereShape>) {
                    computedMass = densityUsed * computeVolumeSphere(s.radius);
                } else if constexpr (std::is_same_v<T, ConeShape>) {
                    computedMass = densityUsed * computeVolumeCone(s.radius, s.height);
                } else if constexpr (std::is_same_v<T, PlaneShape>) {
                    computedMass = (real_t)0;
                } else if constexpr (std::is_same_v<T, MeshShape>) {
                    const MeshAsset& ma = system.assets().getMesh(s.path, s.scale, true);
                    if (ma.volume > (real_t)0) computedMass = densityUsed * ma.volume;
                }
            },
            shape);

        if (computedMass > (real_t)0) massOpt = computedMass;
    }

    const real_t mass = std::max((real_t)0, massOpt.value_or((real_t)0));
    const real_t mu = (props.friction >= (real_t)0) ? props.friction : cfg.friction_default_mu;
    reg.emplace<cardillo::C_Friction>(e, cardillo::C_Friction{mu});

    std::visit(
        [&](auto&& s) {
            using T = std::decay_t<decltype(s)>;

            if constexpr (std::is_same_v<T, CubeShape>) {
                if (props.visual) reg.emplace<cardillo::C_CubeVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) {
                    reg.emplace<cardillo::C_RB_Cube>(
                        e, cardillo::C_RB_Cube{Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity()});
                }
                reg.emplace<cardillo::C_Cube>(e, cardillo::C_Cube{Vector3r::Zero(), s.halfExtents, Quaternion4r::Identity()});
                if (mass > (real_t)0) {
                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{boxInertiaDiag(mass, s.halfExtents)});
                }
            } else if constexpr (std::is_same_v<T, SphereShape>) {
                if (props.visual) reg.emplace<cardillo::C_PointVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) reg.emplace<cardillo::C_RB_Sphere>(e);
                reg.emplace<cardillo::C_Radius>(e, cardillo::C_Radius{s.radius});
                if (mass > (real_t)0) {
                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{sphereInertiaDiag(mass, s.radius)});
                }
            } else if constexpr (std::is_same_v<T, ConeShape>) {
                if (props.visual) reg.emplace<cardillo::C_ConeVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) {
                    reg.emplace<cardillo::C_RB_Cone>(e, cardillo::C_RB_Cone{s.radius, s.height});
                }
                reg.emplace<cardillo::C_Cone>(e, cardillo::C_Cone{s.radius, s.height});
                if (mass > (real_t)0) {
                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{coneInertiaDiag(mass, s.radius, s.height)});
                }
            } else if constexpr (std::is_same_v<T, CapsuleShape>) {
                if (props.visual) reg.emplace<cardillo::C_CapsuleVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) {
                    reg.emplace<cardillo::C_RB_Capsule>(e, cardillo::C_RB_Capsule{s.radius, s.halfLength});
                }
                reg.emplace<cardillo::C_Capsule>(e, cardillo::C_Capsule{s.radius, s.halfLength});
                if (mass > (real_t)0) {
                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{capsuleInertiaDiag(mass, s.radius, s.halfLength)});
                }
            } else if constexpr (std::is_same_v<T, CylinderShape>) {
                if (props.visual) reg.emplace<cardillo::C_CylinderVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) {
                    reg.emplace<cardillo::C_RB_Cylinder>(e, cardillo::C_RB_Cylinder{s.radius, s.halfLength});
                }
                reg.emplace<cardillo::C_Cylinder>(e, cardillo::C_Cylinder{s.radius, s.halfLength});
                if (mass > (real_t)0) {
                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    reg.emplace<cardillo::C_InertiaDiag>(
                        e, cardillo::C_InertiaDiag{cylinderInertiaDiag(mass, s.radius, s.halfLength)});
                }
            } else if constexpr (std::is_same_v<T, PlaneShape>) {
                if (props.visual) reg.emplace<cardillo::C_PlaneVisualTag>(e);
                if (props.collidable && !cfg.collision_disable_all) {
                    reg.emplace<cardillo::C_RB_Plane>(
                        e, cardillo::C_RB_Plane{s.normal, s.up, s.sizeX, s.sizeY});
                }
                reg.emplace<cardillo::C_Plane>(e, cardillo::C_Plane{s.normal, s.up, s.sizeX, s.sizeY});
            } else if constexpr (std::is_same_v<T, MeshShape>) {
                if (props.visual) reg.emplace<cardillo::C_MeshVisualTag>(e);
                reg.emplace<cardillo::C_Mesh>(e, cardillo::C_Mesh{s.path, s.scale});

                const bool dynamic = mass > (real_t)0;
                const MeshAsset& asset = system.assets().getMesh(s.path, s.scale, dynamic);

                if (props.collidable && !cfg.collision_disable_all) {
                    if (s.use_bbox_collider) {
                        Vector3r center = Vector3r::Zero();
                        Vector3r he = Vector3r((real_t)0.05, (real_t)0.05, (real_t)0.05);
                        meshAabbFromAsset(asset, center, he);
                        Quaternion4r q_local(asset.Rpa.transpose());
                        reg.emplace<cardillo::C_RB_Cube>(
                            e, cardillo::C_RB_Cube{center, he, q_local});
                        if (s.show_collider) {
                            reg.emplace<cardillo::C_Cube>(e, cardillo::C_Cube{center, he, q_local});
                            reg.emplace<cardillo::C_CubeVisualTag>(e);
                        }
                    } else {
                        reg.emplace<cardillo::C_RB_Mesh>(e);
                    }
                }

                if (dynamic) {
                    const Quaternion4r q_rpa(asset.Rpa);
                    const Quaternion4r q_new = state.orientation * q_rpa;
                    const Vector3r pos_new = state.position + (state.orientation * asset.com);

                    reg.get<cardillo::C_Position3>(e).value = pos_new;
                    reg.get<cardillo::C_Orientation>(e).value = Quaternion4r(q_new).normalized();

                    reg.emplace<cardillo::C_PhysicsObject>(e);
                    reg.emplace<cardillo::C_RigidBodyTag>(e);
                    reg.emplace<cardillo::C_Mass>(e, cardillo::C_Mass{mass});
                    if (asset.volume > (real_t)0) {
                        const real_t rho = mass / asset.volume;
                        const Vector3r Idiag = rho * asset.inertia_diag_unit.cwiseMax(Vector3r::Zero());
                        reg.emplace<cardillo::C_InertiaDiag>(e, cardillo::C_InertiaDiag{Idiag});
                    }
                }
            }
        },
        shape);

    system.markStructureDirty();
    return e;
}

} // namespace cardillo::physics
