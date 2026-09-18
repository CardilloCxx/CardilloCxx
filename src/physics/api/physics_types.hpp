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
    real_t width{0};   // for rectangular cross-section
    real_t height{0};  // for rectangular cross-section
    real_t radius{0};  // for round cross-section
    std::vector<Vector2r> polygon;  // for polygon cross-section, always centered on its own centroid
    BeamCrossSectionType type{BeamCrossSectionType::Square};

    // Computes the centroid of a simple polygon (convex or concave) via Green's theorem / shoelace formula.
    static Vector2r centroidOf(const std::vector<Vector2r>& poly) {
        real_t A = 0, Cx = 0, Cy = 0;
        size_t n = poly.size();
        for (size_t i = 0; i < n; ++i) {
            const Vector2r& p1 = poly[i];
            const Vector2r& p2 = poly[(i + 1) % n];
            real_t cross = p1.x() * p2.y() - p2.x() * p1.y();
            A += cross;
            Cx += (p1.x() + p2.x()) * cross;
            Cy += (p1.y() + p2.y()) * cross;
        }
        A *= (real_t)0.5;
        if (std::abs(A) < (real_t)1e-12) {
            // Degenerate polygon (zero area) — fall back to vertex average to avoid div-by-zero.
            Vector2r avg(0, 0);
            for (const auto& p : poly) avg += p;
            return n > 0 ? avg / (real_t)n : avg;
        }
        return Vector2r(Cx / ((real_t)6.0 * A), Cy / ((real_t)6.0 * A));
    }

    // Shifts every vertex so the polygon's centroid sits at the origin.
    static std::vector<Vector2r> recenter(const std::vector<Vector2r>& poly) {
        Vector2r c = centroidOf(poly);
        std::vector<Vector2r> out;
        out.reserve(poly.size());
        for (const auto& p : poly) out.push_back(p - c);
        return out;
    }

    static BeamCrossSection square(real_t w, real_t h) {
        BeamCrossSection sec;
        sec.width = w;
        sec.height = h;
        sec.type = BeamCrossSectionType::Square;
        sec.polygon = recenter({Vector2r(-w / 2, -h / 2), Vector2r(w / 2, -h / 2),
                                 Vector2r(w / 2, h / 2), Vector2r(-w / 2, h / 2)});
        return sec;
    }

    static BeamCrossSection triangle(real_t w, real_t h) {
        BeamCrossSection sec;
        sec.width = w;
        sec.height = h;
        sec.type = BeamCrossSectionType::Triangle;
        // NOTE: recentered on centroid, so this is no longer the same local frame as the
        // old vertex layout (base at y=-h/2, apex at y=+h/2) — the *shape* is identical,
        // just shifted so the centroid (not the base) sits at the origin.
        sec.polygon = recenter({Vector2r(-w / 2, -h / 2), Vector2r(w / 2, -h / 2), Vector2r(0, h / 2)});
        return sec;
    }

    static BeamCrossSection round(real_t radius, size_t numSegments = 8) {
        BeamCrossSection sec;
        sec.radius = radius;
        sec.type = BeamCrossSectionType::Round;
        sec.polygon.reserve(numSegments);
        for (size_t i = 0; i < numSegments; ++i) {
            real_t angle = (real_t)i / (real_t)numSegments * (real_t)2.0 * (real_t)M_PI;
            sec.polygon.push_back(Vector2r(radius * std::cos(angle), radius * std::sin(angle)));
        }
        // Already centered by symmetry, but recenter anyway for consistency/robustness
        // (e.g. low numSegments or future changes to vertex generation).
        sec.polygon = recenter(sec.polygon);
        return sec;
    }

    static BeamCrossSection custom(const std::vector<Vector2r>& poly) {
        BeamCrossSection sec;
        sec.polygon = recenter(poly);
        sec.type = BeamCrossSectionType::Polygon;
        return sec;
    }

    real_t area() const {
        switch (type) {
            case BeamCrossSectionType::Round:
                return (real_t)M_PI * radius * radius;
            case BeamCrossSectionType::Triangle:
                return (real_t)0.5 * width * height;
            case BeamCrossSectionType::Polygon: {
                real_t a = (real_t)0;
                size_t n = polygon.size();
                for (size_t i = 0; i < n; ++i) {
                    const Vector2r& p1 = polygon[i];
                    const Vector2r& p2 = polygon[(i + 1) % n];
                    a += p1.x() * p2.y() - p2.x() * p1.y();
                }
                return std::abs(a) * (real_t)0.5;
            }
            default:
                return width * height;
        }
    }

    // Exact second moment of area about the (centroid-aligned) y-axis, for any simple polygon.
    real_t Iy() const {
        switch (type) {
            case BeamCrossSectionType::Round:
                return (real_t)M_PI * std::pow(radius, 4) / (real_t)4.0;
            case BeamCrossSectionType::Triangle:
                return width * std::pow(height, (real_t)3) / (real_t)36.0;
            case BeamCrossSectionType::Polygon: {
                real_t Iy = (real_t)0;
                size_t n = polygon.size();
                for (size_t i = 0; i < n; ++i) {
                    const Vector2r& p1 = polygon[i];
                    const Vector2r& p2 = polygon[(i + 1) % n];
                    Iy += (p1.x() * p2.y() - p2.x() * p1.y()) * (p1.y() * p1.y() + p1.y() * p2.y() + p2.y() * p2.y());
                }
                return std::abs(Iy) / (real_t)12.0;
            }
            default:
                return width * std::pow(height, (real_t)3) / (real_t)12.0;
        }
    }

    // Exact second moment of area about the (centroid-aligned) z-axis, for any simple polygon.
    real_t Iz() const {
        switch (type) {
            case BeamCrossSectionType::Round:
                return (real_t)M_PI * std::pow(radius, 4) / (real_t)4.0;
            case BeamCrossSectionType::Triangle:
                return std::pow(width, (real_t)3) * height / (real_t)36.0;
            case BeamCrossSectionType::Polygon: {
                real_t Iz = (real_t)0;
                size_t n = polygon.size();
                for (size_t i = 0; i < n; ++i) {
                    const Vector2r& p1 = polygon[i];
                    const Vector2r& p2 = polygon[(i + 1) % n];
                    Iz += (p1.x() * p2.y() - p2.x() * p1.y()) * (p1.x() * p1.x() + p1.x() * p2.x() + p2.x() * p2.x());
                }
                return std::abs(Iz) / (real_t)12.0;
            }
            default:
                return std::pow(width, (real_t)3) * height / (real_t)12.0;
        }
    }

    real_t Jp() const { return Iy() + Iz(); }

    // Max |y| / |x| among vertices — the governing extreme-fiber distance from the
    // centroid, used for a single (conservative) section modulus value.
    real_t maxAbsY() const {
        real_t m = (real_t)0;
        for (const auto& p : polygon) m = std::max(m, std::abs(p.y()));
        return m;
    }

    real_t maxAbsX() const {
        real_t m = (real_t)0;
        for (const auto& p : polygon) m = std::max(m, std::abs(p.x()));
        return m;
    }

    real_t sectionModulus() const {
        switch (type) {
            case BeamCrossSectionType::Round:
                return (real_t)M_PI * std::pow(radius, 3) / (real_t)4.0;
            case BeamCrossSectionType::Triangle: {
                real_t S_horizontal = width * height * height / (real_t)24.0;   // apex-side, governs vertical bending
                real_t S_vertical   = width * width * height / (real_t)18.0;    // symmetric, horizontal bending
                return std::min(S_horizontal, S_vertical);
            }
            case BeamCrossSectionType::Polygon: {
                real_t cy = maxAbsY();
                real_t cz = maxAbsX();
                real_t Wy = (cy > (real_t)0) ? Iy() / cy : (real_t)0;
                real_t Wz = (cz > (real_t)0) ? Iz() / cz : (real_t)0;
                return std::min(Wy, Wz);
            }
            default:
                return std::max(width, height) * std::min(width, height) * std::min(width, height) / (real_t)6.0;
        }
    }
};

/// A beam hull shape, used for collision detection between beam segments.
struct BeamHullShape {
    BeamCrossSection cross_section;
    float length{0};

    BeamHullShape() = default;
    explicit BeamHullShape(const BeamCrossSection& cs, float len) : cross_section(cs),length(len) {}
};

using RigidShape = std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape, BeamHullShape>;

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
