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
    /// World-frame linear velocity (m/s).
    Vector3r linearVelocity = Vector3r::Zero();
    /// Body-frame angular velocity (rad/s).
    Vector3r angularVelocity = Vector3r::Zero();

    static RigidState inertial() { return RigidState{}; }

    RigidState() = default;
    explicit RigidState(const Vector3r& p) : position(p) {}
    RigidState(const Vector3r& p, const Vector3r& v) : position(p), linearVelocity(v) {}
    RigidState(const Vector3r& p, const Vector3r& v, const Vector3r& w)
        : position(p), linearVelocity(v), angularVelocity(w) {}
    RigidState(const Vector3r& p, const Quaternion4r& q) : position(p) { setOrientation(q); }
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q)
        : position(p), linearVelocity(v) { setOrientation(q); }
    RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q, const Vector3r& w)
        : position(p), linearVelocity(v), angularVelocity(w) { setOrientation(q); }

    RigidState(const Vector3r& p_local, const Vector3r& v_local, const Quaternion4r& q_local,
               const Vector3r& w_local, entt::entity refEntity, entt::registry& reg);
    RigidState(const Vector3r& p_local, entt::entity refEntity, entt::registry& reg)
        : RigidState(p_local, Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}
    RigidState(entt::entity refEntity, entt::registry& reg)
        : RigidState(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}

    /// World-frame orientation (always unit-length).
    const Quaternion4r& orientation() const { return orientation_; }
    /// Rotation matrix, always consistent with orientation().
    const Matrix33r& rotation() const { return rotation_; }

    void setOrientation(const Quaternion4r& q) {
        Quaternion4r qn = q;
        if (!qn.coeffs().allFinite() ||
            qn.coeffs().squaredNorm() <= std::numeric_limits<real_t>::epsilon()) {
            qn = Quaternion4r::Identity();
        } else {
            qn.normalize();
        }
        orientation_ = qn;
        rotation_ = qn.toRotationMatrix();
    }

    void setRotation(const Matrix33r& R) { setOrientation(Quaternion4r(R)); }

private:
    Quaternion4r orientation_ = Quaternion4r::Identity();
    Matrix33r rotation_ = Matrix33r::Identity();
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
    Matrix33r localRotation{Matrix33r::Identity()};
    CapsuleShape() = default;
    CapsuleShape(real_t r, real_t h, const Matrix33r& localR = Matrix33r::Identity()) : radius(r), halfLength(h), localRotation(localR) {}
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

enum class BeamColliderMode { RigidBodyPrimitive, InterSegmentHull };
enum class BeamCrossSectionType { Square, Triangle, Round, Polygon };

struct BeamCrossSection {
    real_t width{0};   
    real_t height{0};  
    real_t radius{0};  
    std::vector<Vector2r> polygon;  
    BeamCrossSectionType type{BeamCrossSectionType::Square};

    static Vector2r centroidOf(const std::vector<Vector2r>& poly);
    static std::vector<Vector2r> recenter(const std::vector<Vector2r>& poly);

    static BeamCrossSection rectangular(real_t w, real_t h);
    static BeamCrossSection triangle(real_t w, real_t h);
    static BeamCrossSection round(real_t radius, size_t numSegments = 8);
    static BeamCrossSection custom(const std::vector<Vector2r>& poly);

    real_t area() const;
    real_t Iy() const;
    real_t Iz() const;
    real_t Jp() const;
    real_t maxAbsY() const;
    real_t maxAbsX() const;
    real_t sectionModulus() const;
};

struct BeamHullShape {
    BeamCrossSection cross_section;
    float length{0};

    BeamHullShape() = default;
    explicit BeamHullShape(const BeamCrossSection& cs, float len);
};

using RigidShape = std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape, BeamHullShape>;

struct BeamSpringParams {
    real_t E{0};
    real_t nu{0};
    Vector3r scaleKe{Vector3r::Ones()};
    Vector3r scaleKf{Vector3r::Ones()};
    std::optional<Vector3r> Ke_direct;
    std::optional<Vector3r> Kf_direct;
    std::optional<Vector3r> gamma0;
    std::optional<Vector3r> kappa0;
    real_t dampingFactor = 0.0;

    BeamSpringParams() = default;
    BeamSpringParams(const Vector3r& Ke_in, const Vector3r& Kf_in, real_t dampingFactor_in = 0.0);

    Vector3r Ke(real_t segLen, const BeamCrossSection& sec) const;
    Vector3r Kf(real_t segLen, const BeamCrossSection& sec) const;
    void setDampingFromFactor(real_t d);

    static BeamSpringParams fromMaterial(real_t E_in, real_t nu_in, real_t axialScale = (real_t)1, 
                                          real_t shearScale = (real_t)1, real_t torsionScale = (real_t)1, 
                                          real_t bendYScale = (real_t)1, real_t bendZScale = (real_t)1, 
                                          real_t dampingFactor_in = (real_t)0);
};

}  // namespace cardillo::physics
