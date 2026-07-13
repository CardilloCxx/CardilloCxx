#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <variant>

#include <entt/entt.hpp>

#include "../../misc/types.hpp"

namespace cardillo {
class World;
}

namespace cardillo::physics {

/// Initial kinematic state for a rigid body: pose and velocities in the world frame.
struct RigidState {
    /// World-frame centre-of-mass position (metres).
    Vector3r position = Vector3r::Zero();
    /// World-frame orientation as a unit quaternion.
    Quaternion4r orientation = Quaternion4r::Identity();
    /// Rotation matrix derived from @p orientation. The named constructors keep this in sync;
    /// if you set @p orientation directly you must update @p rotation as well.
    Matrix33r rotation = Matrix33r::Identity();
    /// World-frame linear velocity (m/s).
    Vector3r linearVelocity = Vector3r::Zero();
    /// Body-frame angular velocity (rad/s).
    Vector3r angularVelocity = Vector3r::Zero();

    static RigidState inertial() { return RigidState{}; }

    RigidState() = default;
    explicit RigidState(const Vector3r& p) : position(p) {}
    RigidState(const Vector3r& p, const Vector3r& v) : position(p), linearVelocity(v) {}
    RigidState(const Vector3r& p, const Quaternion4r& q) : position(p), orientation(q), rotation(q.toRotationMatrix()) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q) : position(p), orientation(q), rotation(q.toRotationMatrix()), linearVelocity(v) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q, const Vector3r& w) : position(p), orientation(q), rotation(q.toRotationMatrix()), linearVelocity(v), angularVelocity(w) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Vector3r& w) : position(p), linearVelocity(v), angularVelocity(w) {}

    RigidState(const Vector3r& p_local, const Vector3r& v_local, const Quaternion4r& q_local, const Vector3r& w_local, entt::entity refEntity, entt::registry& reg);

    RigidState(const Vector3r& p_local, entt::entity refEntity, entt::registry& reg) : RigidState(p_local, Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}

    RigidState(entt::entity refEntity, entt::registry& reg) : RigidState(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}
};

/// Axis-aligned box shape. Half-extents are in the body's local frame.
struct CubeShape {
    /// Half-extents along local x, y, z axes (metres).
    Vector3r halfExtents{Vector3r::Zero()};
    CubeShape() = default;
    explicit CubeShape(const Vector3r& he) : halfExtents(he) {}
};

/// Infinite flat surface (collision-wise). The visual quad is finite.
/// Use with @ref PhysicsEngine::addStaticBody — PlaneShape has no inertia.
struct PlaneShape {
    /// Outward surface normal in world space.
    Vector3r normal{Vector3r(0, 0, 1)};
    /// Up direction for the visual quad.
    Vector3r up{Vector3r(0, 1, 0)};
    /// Visual half-width along the tangent axis (metres).
    real_t sizeX{5};
    /// Visual half-width along the up axis (metres).
    real_t sizeY{5};
    PlaneShape() = default;
    PlaneShape(const Vector3r& n, const Vector3r& u, real_t sx, real_t sy) : normal(n), up(u), sizeX(sx), sizeY(sy) {}
};

/// Cylinder capped with two hemispheres. Long axis runs along the body's local z-axis.
struct CapsuleShape {
    /// Hemisphere cap radius (metres). Total extent along local z: ±(halfLength + radius).
    real_t radius{0};
    /// Half-length of the cylindrical shaft between the two caps (metres).
    real_t halfLength{0};
    CapsuleShape() = default;
    CapsuleShape(real_t r, real_t h) : radius(r), halfLength(h) {}
};

/// Flat-ended cylinder. Long axis runs along the body's local z-axis.
struct CylinderShape {
    /// Barrel radius (metres).
    real_t radius{0};
    /// Half of the total cylinder height (metres).
    real_t halfLength{0};
    CylinderShape() = default;
    CylinderShape(real_t r, real_t h) : radius(r), halfLength(h) {}
};

/// Right circular cone. Tip points along the body's local +z axis.
struct ConeShape {
    /// Base radius (metres).
    real_t radius{0};
    /// Full height from base to tip (metres).
    real_t height{0};
    ConeShape() = default;
    ConeShape(real_t r, real_t h) : radius(r), height(h) {}
};

/// Sphere shape.
struct SphereShape {
    /// Sphere radius (metres).
    real_t radius{0};
    SphereShape() = default;
    explicit SphereShape(real_t r) : radius(r) {}
};

/// Triangle mesh shape loaded from an OBJ or STL file.
/// The mesh is normalised to its principal-axes frame and volume-weighted CoM.
struct MeshShape {
    /// Path to the OBJ or STL file on disk.
    std::string path;
    /// Per-axis scale applied to the mesh on load. Use (0.001,0.001,0.001) for mm→m conversion.
    Vector3r scale{1, 1, 1};
    /// When true, replace the exact mesh hull with its axis-aligned bounding box for collision.
    bool use_bbox_collider{false};
    /// When true, also render the bounding box in the VTK output (only meaningful with use_bbox_collider).
    bool show_collider{false};
    MeshShape() = default;
    explicit MeshShape(const std::string& p, bool bbox = false, bool showCol = false) : path(p), use_bbox_collider(bbox), show_collider(showCol) {}
    MeshShape(const std::string& p, const Vector3r& s, bool bbox = false, bool showCol = false) : path(p), scale(s), use_bbox_collider(bbox), show_collider(showCol) {}
};

using RigidShape = std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape>;

/// Physical properties and pipeline flags for a rigid body.
struct RigidProps {
    /// Body mass (kg). Takes priority over @p density when both are set.
    /// If neither is set the body is created with zero mass and treated as static.
    std::optional<real_t> mass;
    /// Density (kg/m³). Multiplied by the shape volume to compute mass when @p mass is unset.
    std::optional<real_t> density;
    /// Coulomb friction coefficient. Negative value (-1) means use Config::friction_default_mu.
    real_t friction = -1;
    /// Normal coefficient of restitution. Negative value (-1) means use Config::restitution_default_normal.
    real_t restitution_normal = -1;
    /// Tangential coefficient of restitution. Negative value (-1) means use Config::restitution_default_tangential.
    real_t restitution_tangential = -1;
    /// Register this body with the collision detection system.
    bool collidable = true;
    /// Include this body in VTK output.
    bool visual = true;

    RigidProps() = default;
    explicit RigidProps(real_t m) : mass(m) {}
    RigidProps(real_t m, real_t fric, bool vis = true, bool coll = true) : mass(m), friction(fric), collidable(coll), visual(vis) {}

    static RigidProps withDensity(real_t rho) {
        RigidProps p;
        p.density = rho;
        return p;
    }
};

enum class BeamBodyType { Cube, Capsule, Cylinder };

/// Cross-section geometry of a beam segment.
/// Provides derived section properties (area, moments of inertia) used to compute beam stiffness.
struct BeamCrossSection {
    /// Cross-section width (metres). For Capsule/Cylinder types, treated as diameter.
    real_t width{0};
    /// Cross-section height (metres). For Capsule/Cylinder types, treated as diameter.
    real_t height{0};
    /// Shape type used for each beam segment's collision/visual geometry.
    BeamBodyType type{BeamBodyType::Cube};

    BeamCrossSection() = default;
    BeamCrossSection(real_t w, real_t h, BeamBodyType t = BeamBodyType::Cube) : width(w), height(h), type(t) {}

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

/// Elastic and damping parameters for a Cosserat-rod beam constraint.
/// Stiffness can be derived from material constants (E, nu) or set directly via Ke_direct/Kf_direct.
struct BeamSpringParams {
    /// Young's modulus (Pa). Used with @p nu to compute Ke/Kf from cross-section geometry.
    real_t E{0};
    /// Poisson's ratio. Used to compute shear modulus G = E / (2*(1+nu)).
    real_t nu{0};
    /// Per-component scale factors for the extensional/shear stiffness Ke = (E*A/L, G*A/L, G*A/L).
    /// Components: [axial, shear-y, shear-z].
    Vector3r scaleKe{Vector3r::Ones()};
    /// Per-component scale factors for the torsion/bending stiffness Kf = (G*Jp/L, E*Iy/L, E*Iz/L).
    /// Components: [torsion, bend-y, bend-z].
    Vector3r scaleKf{Vector3r::Ones()};
    /// Direct override for Ke (N/m). When set, E, nu, and scaleKe are ignored for stretch/shear.
    std::optional<Vector3r> Ke_direct;
    /// Direct override for Kf (Nm/rad). When set, E, nu, and scaleKf are ignored for torsion/bending.
    std::optional<Vector3r> Kf_direct;
    /// Optional rest-state translational strain γ₀. Initialized from the initial pose when unset.
    std::optional<Vector3r> gamma0;
    /// Optional rest-state curvature κ₀. Initialized from the initial relative orientation when unset.
    std::optional<Vector3r> kappa0;
    /// Rayleigh-type damping factor d. Effective damping stiffness = K * d (units: s).
    real_t dampingFactor = 0.0;

    BeamSpringParams() = default;

    BeamSpringParams(const Vector3r& Ke_in, const Vector3r& Kf_in, real_t dampingFactor_in = 0.0) : Ke_direct(Ke_in), Kf_direct(Kf_in), dampingFactor(dampingFactor_in) {}

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

    static BeamSpringParams fromMaterial(real_t E_in, real_t nu_in, real_t axialScale = (real_t)1, real_t shearScale = (real_t)1, real_t torsionScale = (real_t)1, real_t bendYScale = (real_t)1,
                                         real_t bendZScale = (real_t)1, real_t dampingFactor_in = (real_t)0) {
        BeamSpringParams p;
        p.E = E_in;
        p.nu = nu_in;
        p.scaleKe = Vector3r(axialScale, shearScale, shearScale);
        p.scaleKf = Vector3r(torsionScale, bendYScale, bendZScale);
        p.dampingFactor = dampingFactor_in;
        return p;
    }
};

}  // namespace cardillo::physics
