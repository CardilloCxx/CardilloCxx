#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <variant>

#include <entt/entt.hpp>

#include "../misc/types.hpp"

namespace cardillo {
class World;
}

namespace cardillo::physics {

struct RigidState {
    Vector3r position = Vector3r::Zero();
    Quaternion4r orientation = Quaternion4r::Identity();
    Vector3r linearVelocity = Vector3r::Zero();
    Vector3r angularVelocity = Vector3r::Zero();

    RigidState() = default;
    explicit RigidState(const Vector3r& p) : position(p) {}
    RigidState(const Vector3r& p, const Vector3r& v) : position(p), linearVelocity(v) {}
    RigidState(const Vector3r& p, const Quaternion4r& q) : position(p), orientation(q) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q)
        : position(p), orientation(q), linearVelocity(v) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q, const Vector3r& w)
        : position(p), orientation(q), linearVelocity(v), angularVelocity(w) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Vector3r& w)
        : position(p), linearVelocity(v), angularVelocity(w) {}

    RigidState(const Vector3r& p_local,
               const Vector3r& v_local,
               const Quaternion4r& q_local,
               const Vector3r& w_local,
               entt::entity refEntity,
               entt::registry& reg);

    RigidState(const Vector3r& p_local,
               entt::entity refEntity,
               entt::registry& reg)
        : RigidState(p_local,
                     Vector3r::Zero(),
                     Quaternion4r::Identity(),
                     Vector3r::Zero(),
                     refEntity,
                     reg) {}

    RigidState(entt::entity refEntity,
               entt::registry& reg)
        : RigidState(Vector3r::Zero(),
                     Vector3r::Zero(),
                     Quaternion4r::Identity(),
                     Vector3r::Zero(),
                     refEntity,
                     reg) {}
};

struct CubeShape {
    Vector3r halfExtents{Vector3r::Zero()};
    CubeShape() = default;
    explicit CubeShape(const Vector3r& he) : halfExtents(he) {}
};

struct PlaneShape {
    Vector3r normal{Vector3r(0, 0, 1)};
    Vector3r up{Vector3r(0, 1, 0)};
    real_t sizeX{5};
    real_t sizeY{5};
    PlaneShape() = default;
    PlaneShape(const Vector3r& n, const Vector3r& u, real_t sx, real_t sy)
        : normal(n), up(u), sizeX(sx), sizeY(sy) {}
};

struct CapsuleShape {
    real_t radius{0};
    real_t halfLength{0};
    CapsuleShape() = default;
    CapsuleShape(real_t r, real_t h) : radius(r), halfLength(h) {}
};

struct CylinderShape {
    real_t radius{0};
    real_t halfLength{0};
    CylinderShape() = default;
    CylinderShape(real_t r, real_t h) : radius(r), halfLength(h) {}
};

struct ConeShape {
    real_t radius{0};
    real_t height{0};
    ConeShape() = default;
    ConeShape(real_t r, real_t h) : radius(r), height(h) {}
};

struct SphereShape {
    real_t radius{0};
    SphereShape() = default;
    explicit SphereShape(real_t r) : radius(r) {}
};

struct MeshShape {
    std::string path;
    Vector3r scale{1, 1, 1};
    bool use_bbox_collider{false};
    bool show_collider{false};
    MeshShape() = default;
    explicit MeshShape(const std::string& p, bool bbox = false, bool showCol = false)
        : path(p), use_bbox_collider(bbox), show_collider(showCol) {}
    MeshShape(const std::string& p, const Vector3r& s, bool bbox = false, bool showCol = false)
        : path(p), scale(s), use_bbox_collider(bbox), show_collider(showCol) {}
};

using RigidShape = std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape>;

struct RigidProps {
    std::optional<real_t> mass;
    std::optional<real_t> density;
    real_t friction = -1;
    bool collidable = true;
    bool visual = true;

    RigidProps() = default;
    explicit RigidProps(real_t m) : mass(m) {}
    RigidProps(real_t m, real_t fric, bool vis = true, bool coll = true)
        : mass(m), friction(fric), collidable(coll), visual(vis) {}

    static RigidProps withDensity(real_t rho) {
        RigidProps p;
        p.density = rho;
        return p;
    }
};

enum class BeamBodyType { Cube, Capsule, Cylinder };

struct BeamCrossSection {
    real_t width{0};
    real_t height{0};
    BeamBodyType type{BeamBodyType::Cube};

    BeamCrossSection() = default;
    BeamCrossSection(real_t w, real_t h, BeamBodyType t = BeamBodyType::Cube)
        : width(w), height(h), type(t) {}

    real_t area() const {
        if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
            real_t r = (std::min(width, height)) * (real_t)0.5;
            return (real_t)M_PI * r * r;
        }
        return width * height;
    }

    real_t Iy() const {
        if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
            real_t r = (std::min(width, height)) * (real_t)0.5;
            return (real_t)M_PI * std::pow(r, 4) / (real_t)4.0;
        }
        return width * std::pow(height, (real_t)3) / (real_t)12.0;
    }

    real_t Iz() const {
        if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
            real_t r = (std::min(width, height)) * (real_t)0.5;
            return (real_t)M_PI * std::pow(r, 4) / (real_t)4.0;
        }
        return std::pow(width, (real_t)3) * height / (real_t)12.0;
    }

    real_t Jp() const { return Iy() + Iz(); }

    real_t sectionModulus() const {
        if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
            real_t r = (std::min(width, height)) * (real_t)0.5;
            if (r == (real_t)0.0) return (real_t)0.0;
            return (real_t)M_PI * std::pow(r, 3) / (real_t)4.0;
        }
        real_t Wy = Iy() / ((real_t)0.5 * height);
        real_t Wz = Iz() / ((real_t)0.5 * width);
        return std::min(Wy, Wz);
    }
};

struct BeamSpringParams {
    real_t E{0};
    real_t nu{0};
    Vector3r scaleKe{Vector3r::Ones()};
    Vector3r scaleKf{Vector3r::Ones()};
    std::optional<Vector3r> Ke_direct;
    std::optional<Vector3r> Kf_direct;
    real_t dampingFactor = 0.0;

    BeamSpringParams() = default;

    BeamSpringParams(const Vector3r& Ke_in,
                     const Vector3r& Kf_in,
                     real_t dampingFactor_in = 0.0)
        : Ke_direct(Ke_in), Kf_direct(Kf_in), dampingFactor(dampingFactor_in) {}

    Vector3r Ke(real_t segLen, const BeamCrossSection& sec) const {
        if (Ke_direct.has_value()) return *Ke_direct;
        const real_t G = E / ((real_t)2.0 * ((real_t)1.0 + nu));
        const real_t A = sec.area();
        Vector3r base(E * A / segLen, G * A / segLen, G * A / segLen);
        return base.cwiseProduct(scaleKe);
    }

    Vector3r Kf(real_t segLen, const BeamCrossSection& sec) const {
        if (Kf_direct.has_value()) return *Kf_direct;
        const real_t G = E / ((real_t)2.0 * ((real_t)1.0 + nu));
        Vector3r base(G * sec.Jp() / segLen, E * sec.Iy() / segLen, E * sec.Iz() / segLen);
        return base.cwiseProduct(scaleKf);
    }

    void setDampingFromFactor(real_t d) { dampingFactor = d; }

    static BeamSpringParams fromMaterial(real_t E_in,
                                         real_t nu_in,
                                         real_t axialScale = (real_t)1,
                                         real_t shearScale = (real_t)1,
                                         real_t torsionScale = (real_t)1,
                                         real_t bendYScale = (real_t)1,
                                         real_t bendZScale = (real_t)1,
                                         real_t dampingFactor_in = (real_t)0) {
        BeamSpringParams p;
        p.E = E_in;
        p.nu = nu_in;
        p.scaleKe = Vector3r(axialScale, shearScale, shearScale);
        p.scaleKf = Vector3r(torsionScale, bendYScale, bendZScale);
        p.dampingFactor = dampingFactor_in;
        return p;
    }
};

} // namespace cardillo::physics
