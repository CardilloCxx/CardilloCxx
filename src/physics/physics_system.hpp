#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <optional>
#include <variant>
#include <functional>
#include <unordered_map>
#include <type_traits>
#include <entt/entt.hpp>
#include <petscsys.h>
#include "../misc/types.hpp"
#include "../misc/spline.hpp"
#include "assets.hpp"
#include "../misc/dofs.hpp"
#include "../config/config.hpp"
#include "../misc/timings/TimingManager.hpp"
// COAL types for mesh assets
// match installed COAL include layout (lowercase paths)
#include <coal/BVH/BVH_model.h>
#include <coal/collision_object.h>
#include <coal/broadphase/broadphase.h>
// HeightField collider
#include <coal/hfield.h>

// forward-declare warmstart provider interface to avoid include cycles
namespace cardillo { namespace solver { class WarmstartProvider; } }

// fwd
namespace cardillo { namespace collision { class CollisionCoal; } }

// forward-declare new constraint pattern types (defined in constraints.hpp)
namespace cardillo { namespace physics { class ConstraintPattern; class LinearDistanceConstraint; } }

namespace cardillo {

// A minimal, standard-C++ physics system for frictionless point masses
// with translational DOFs only (no rotations).
class PhysicsSystem {
public:
    static Quaternion4r alignQuaternionTo(const Quaternion4r& q_in, const Quaternion4r& q_ref) {
        Quaternion4r q = q_in;
        if (!q.coeffs().allFinite()) return q_ref;
        q.normalize();
        if (q_ref.dot(q) < (real_t)0) q.coeffs() = -q.coeffs();
        return q;
    }

    // Unified rigid-body creation API --------------------------------------------------
    struct RigidState {
        Vector3r     position       = Vector3r::Zero();
        Quaternion4r orientation    = Quaternion4r::Identity();
        Vector3r     linearVelocity = Vector3r::Zero();
        Vector3r     angularVelocity= Vector3r::Zero();
        RigidState() = default;
        explicit RigidState(const Vector3r& p) : position(p) {}
        RigidState(const Vector3r& p, const Vector3r& v) : position(p), linearVelocity(v) {}
        RigidState(const Vector3r& p, const Quaternion4r& q) : position(p), orientation(q) {}
        RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q)
            : position(p), orientation(q), linearVelocity(v) {}
        RigidState(const Vector3r& p, const Vector3r& v, const Quaternion4r& q, const Vector3r& w)
            : position(p), orientation(q), linearVelocity(v), angularVelocity(w) {}
        // Convenience when no rotation is desired but angular velocity is provided
        RigidState(const Vector3r& p, const Vector3r& v, const Vector3r& w)
            : position(p), linearVelocity(v), angularVelocity(w) {}
        // Convenience: full state specified in a moving reference frame.
        // p_local, v_local, q_local, w_local are expressed in the reference frame.
        // They are transformed into world coordinates using the reference's pose and velocities.
        RigidState(const Vector3r& p_local,
                   const Vector3r& v_local,
                   const Quaternion4r& q_local,
                   const Vector3r& w_local,
                   entt::entity refEntity,
                   entt::registry& reg) {
            if (refEntity != entt::null &&
                reg.all_of<PhysicsSystem::C_Position3,
                           PhysicsSystem::C_Orientation,
                           PhysicsSystem::C_LinearVelocity3,
                           PhysicsSystem::C_AngularVelocity3>(refEntity)) {

                const auto& r_ORef = reg.get<PhysicsSystem::C_Position3>(refEntity).value;      // origin of ref in world
                const auto& q_Ref  = reg.get<PhysicsSystem::C_Orientation>(refEntity).value;     // ref orientation in world
                const auto& v_Ref  = reg.get<PhysicsSystem::C_LinearVelocity3>(refEntity).value; // linear vel of ref in world
                const auto& w_Ref  = reg.get<PhysicsSystem::C_AngularVelocity3>(refEntity).value;// angular vel of ref in world

                const Matrix33r A_Ref = q_Ref.toRotationMatrix();
                position    = r_ORef + A_Ref * p_local;
                orientation = q_Ref * q_local;

                const Vector3r r_rel_world = A_Ref * p_local;
                linearVelocity  = v_Ref + w_Ref.cross(r_rel_world) + A_Ref * v_local;
                angularVelocity = w_Ref + A_Ref * w_local;
            } else {
                position        = p_local;
                orientation     = q_local;
                linearVelocity  = v_local;
                angularVelocity = w_local;
            }
        }

         RigidState(const Vector3r& p_local,
                   entt::entity refEntity,
                   entt::registry& reg) : RigidState(p_local, Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}
        
        RigidState(entt::entity refEntity,
                   entt::registry& reg) : RigidState(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity(), Vector3r::Zero(), refEntity, reg) {}
    };
    struct CubeShape    { 
        Vector3r halfExtents{Vector3r::Zero()};
        CubeShape() = default;
        explicit CubeShape(const Vector3r& he) : halfExtents(he) {}
    };
    struct PlaneShape   { 
        Vector3r normal{Vector3r(0,0,1)}; Vector3r up{Vector3r(0,1,0)}; real_t sizeX{5}; real_t sizeY{5};
        PlaneShape() = default;
        PlaneShape(const Vector3r& n, const Vector3r& u, real_t sx, real_t sy) : normal(n), up(u), sizeX(sx), sizeY(sy) {}
    };
    struct CapsuleShape { 
        real_t radius{0}; real_t halfLength{0};
        CapsuleShape() = default;
        CapsuleShape(real_t r, real_t h) : radius(r), halfLength(h) {}
    };
    struct CylinderShape {
        real_t radius{0}; real_t halfLength{0};
        CylinderShape() = default;
        CylinderShape(real_t r, real_t h) : radius(r), halfLength(h) {}
    };
    struct ConeShape {
        real_t radius{0};
        real_t height{0};
        ConeShape() = default;
        ConeShape(real_t r, real_t h) : radius(r), height(h) {}
    };
    struct SphereShape  { 
        real_t radius{0};
        SphereShape() = default;
        explicit SphereShape(real_t r) : radius(r) {}
    };
    struct MeshShape    { 
        std::string path; Vector3r scale{1,1,1}; bool use_bbox_collider{false}; bool show_collider{false};
        MeshShape() = default;
        explicit MeshShape(const std::string& p, bool bbox=false, bool showCol=false) : path(p), use_bbox_collider(bbox), show_collider(showCol) {}
        MeshShape(const std::string& p, const Vector3r& s, bool bbox=false, bool showCol=false) : path(p), scale(s), use_bbox_collider(bbox), show_collider(showCol) {}
    };
    using RigidShape = std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape>;
    struct RigidProps {
        // If neither mass nor density set => static (no physics object)
        std::optional<real_t> mass;     
        std::optional<real_t> density; 
        real_t friction = -1; // <0 => use default from config
        bool   collidable = true;
        bool   visual     = true;
        RigidProps() = default;
        explicit RigidProps(real_t m) : mass(m) {}
        RigidProps(real_t m, real_t fric, bool vis=true, bool coll=true) : mass(m), friction(fric), collidable(coll), visual(vis) {}
        static RigidProps withDensity(real_t rho) { RigidProps p; p.density = rho; return p; }
    };

    // Beam cross-section and spring params --------------------------------------------------
    enum class BeamBodyType { Cube, Capsule, Cylinder };
    struct BeamCrossSection {
        real_t width{0};
        real_t height{0};
        BeamBodyType type{BeamBodyType::Cube};
        BeamCrossSection() = default;
        BeamCrossSection(real_t w, real_t h, BeamBodyType t=BeamBodyType::Cube) : width(w), height(h), type(t) {}
        // Area
        real_t area() const {
            if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
                real_t r = (std::min(width, height)) * (real_t)0.5;
                return (real_t)M_PI * r * r;
            }
            return width * height;
        }
        // Second moments of area about local Y and Z
        real_t Iy() const {
            if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
                real_t r = (std::min(width, height)) * (real_t)0.5;
                return (real_t)M_PI * std::pow(r, 4) / (real_t)4.0; // circle
            }
            return width * std::pow(height, (real_t)3) / (real_t)12.0;
        }
        real_t Iz() const {
            if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
                real_t r = (std::min(width, height)) * (real_t)0.5;
                return (real_t)M_PI * std::pow(r, 4) / (real_t)4.0; // circle
            }
            return std::pow(width, (real_t)3) * height / (real_t)12.0;
        }
        real_t Jp() const { return Iy() + Iz(); } // polar approx

        real_t sectionModulus() const {
        if (type == BeamBodyType::Capsule || type == BeamBodyType::Cylinder) {
            // Circle: c = R, Iy = Iz = I. W = I/R = (pi*R^4/4) / R = pi*R^3/4
            real_t r = (std::min(width, height)) * (real_t)0.5;
            if (r == (real_t)0.0) return (real_t)0.0;
            return (real_t)M_PI * std::pow(r, 3) / (real_t)4.0;
        } else { // Cube/Rectangle
            // W_y = Iy / (height/2) 
            real_t Wy = Iy() / ((real_t)0.5 * height);
            // W_z = Iz / (width/2)
            real_t Wz = Iz() / ((real_t)0.5 * width);
            return std::min(Wy, Wz);
        }
    }
    };

    struct BeamSpringParams {
        // Material properties for derived stiffness
        real_t E{0};
        real_t nu{0};
        // Independent scales for each component
        Vector3r scaleKe{Vector3r::Ones()}; // [axial, shearY, shearZ]
        Vector3r scaleKf{Vector3r::Ones()}; // [torsion, bendY, bendZ]
        // Rankine cracking model: maximum crack strain threshold for shear retention
        real_t crackStrainMax{std::numeric_limits<real_t>::infinity()};
        real_t tensileStrength{std::numeric_limits<real_t>::infinity()};
        // Optional direct per-segment stiffness overrides (units already [*/L])
        std::optional<Vector3r> Ke_direct;
        std::optional<Vector3r> Kf_direct;
        // Damping compliances
        real_t dampingFactor = 0.0;

        BeamSpringParams() = default;
        // Direct constructor for per-segment stiffness vectors
        BeamSpringParams(const Vector3r& Ke_in, const Vector3r& Kf_in,
                         real_t dampingFactor_in = 0.0)
            : Ke_direct(Ke_in), Kf_direct(Kf_in), dampingFactor(dampingFactor_in) {}

        // Accessors to per-segment stiffness vectors given section and segment length
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

        void setDampingFromFactor(real_t d) {
            dampingFactor = d;
        }

        // Factory for material-based parameters without needing section at construction
        static BeamSpringParams fromMaterial(real_t E_in, real_t nu_in,
                                             real_t axialScale=(real_t)1,
                                             real_t shearScale=(real_t)1,
                                             real_t torsionScale=(real_t)1,
                                             real_t bendYScale=(real_t)1,
                                             real_t bendZScale=(real_t)1,
                                             real_t dampingFactor_in=(real_t)0)
        {
            BeamSpringParams p;
            p.E = E_in; p.nu = nu_in;
            p.scaleKe = Vector3r(axialScale, shearScale, shearScale);
            p.scaleKf = Vector3r(torsionScale, bendYScale, bendZScale);
            p.dampingFactor = dampingFactor_in;
            return p;
        }
    };

    PhysicsSystem();
    explicit PhysicsSystem(const config::Config& cfg);
    ~PhysicsSystem();
    // Global config accessible across subsystems
    void setConfig(const config::Config& cfg) { m_cfg = cfg; setGravity(m_cfg.sim_gravity); }
    const config::Config& config() const { return m_cfg; }

    // Persistent collision manager (COAL) storage and access
    collision::CollisionCoal& collisionManager();
    const collision::CollisionCoal& collisionManager() const;
    // Timings access
    cardillo::misc::TimingManager& timings();
    const cardillo::misc::TimingManager& timings() const;


    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const { return m_gravity; }

    // Visual/collision plane configuration (for future use in collisions)
    struct Plane {
        Vector3r center{0,0,0};
        Vector3r normal{0,0,1};
        Vector3r up{0,1,0};
        real_t sizeX{5}, sizeY{5}; // half extents for visualization
    };
    void updateBeamElementEntity(entt::entity e);
    void updateEntities();

    // Dynamics getters (Cache them inside the entity to avoid recomputation)
    MatrixXXr getMass( entt::entity e ) const;        // Linear Inertia and Angular Inertia
    // Inverse mass diagonal (vector) for the entity's velocity-space dofs.
    // Order matches getVelocity():
    //  - Rigid body: [1/m, 1/m, 1/m, 1/Ixx, 1/Iyy, 1/Izz]
    //  - Point mass: [1/m, 1/m, 1/m]
    VectorXr getMassInverseDiag( entt::entity e ) const;
    // Generalized inertia diagonal getter (body frame). Computes or retrieves best-known diag.
    Vector3r getInertiaDiag(entt::entity e) const;
    VectorXr getPosition( entt::entity e ) const;     // Linear and angular combined
    VectorXr getVelocity( entt::entity e ) const;     // Linear and angular combined
    VectorXr getForce( entt::entity e ) const;        // Linear and angular combined (gravity + external + gyroscopic)
    VectorXr getForceExternal( entt::entity e ) const; // Gravity + external forces/torques only
    VectorXr getForceGyroscopic( entt::entity e ) const; // Gyroscopic torque only (tau = -w x (I*w))
    real_t getKineticEnergy( entt::entity e ) const;

    // Shared asset access (wrappers over PhysicsAssets using entity components)
    const ::cardillo::MeshAsset& getMeshAsset(entt::entity e) const;
    const ::cardillo::HeightFieldAsset& getHeightFieldAsset(entt::entity e) const;

    // Optional inertia component for rigid bodies with arbitrary shapes
    struct C_InertiaDiag { Vector3r I; }; // body-frame diagonal inertia (Ixx,Iyy,Izz)

    // Expose ECS for external querying (read-only)
    const entt::registry& ecs() const { return m_reg; }
    // Mutable ECS access when external systems need to emplace components
    entt::registry& ecs() { return m_reg; }

    // Access to warmstart provider owned by the system (may be nullptr)
    cardillo::solver::WarmstartProvider* warmstartProvider() const { return m_warmstart_provider.get(); }

    // Disable collision between two entities (order-independent). This persists until enabled again.
    void disableCollisionBetween(entt::entity a, entt::entity b);

    // New constraint-pattern API -------------------------------------------
    // Access all constraint patterns (mutable and const)
    std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>>& constraintPatterns() { return m_constraints_new; }
    const std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>>& constraintPatterns() const { return m_constraints_new; }

    // Add a general constraint pattern (ownership transferred). Returns index.
    size_t addConstraint(std::unique_ptr<cardillo::physics::ConstraintPattern> pattern);

    // Convenience: construct and add a constraint pattern in-place.
    // Usage: sys.addConstraint<physics::LinearDistanceConstraint>(reg, a, b, rA, rB, k, d);
    template <typename Pattern, typename... Args>
    size_t addConstraint(Args&&... args) {
        static_assert(std::is_base_of<cardillo::physics::ConstraintPattern, Pattern>::value,
                      "Pattern must derive from ConstraintPattern");
        return addConstraint(std::unique_ptr<cardillo::physics::ConstraintPattern>(
            new Pattern(std::forward<Args>(args)...)));
    }

    // Public ECS component/tag types for queries
    struct C_Mass { real_t m; };
    struct C_Position3 { Vector3r value; };
    struct C_LinearVelocity3 { Vector3r value; };
    struct C_AngularVelocity3 { Vector3r value; };
    struct C_Orientation { Quaternion4r value; };
    struct C_DirectorTriad { Matrix33r value; };
    struct C_PhysicsObject {};
    struct C_PointMassTag {};
    struct C_RigidBodyTag {};
    struct C_RigidBodyDirectorTag {};
    struct C_Plane { Vector3r normal; Vector3r up; real_t sizeX; real_t sizeY; };
    struct C_Cube { Vector3r center{Vector3r::Zero()}; Vector3r halfExtents; Quaternion4r q{Quaternion4r::Identity()}; };
    struct C_Capsule { real_t radius; real_t halfLength; };
    struct C_Cylinder { real_t radius; real_t halfLength; };
    struct C_Cone { real_t radius; real_t height; };
    struct C_Friction { real_t mu; }; // optional friction coefficient per entity (>=0), absent => 0
    struct C_VisualObject {};
    struct C_PointVisualTag {};
    struct C_PlaneVisualTag {};
    struct C_CubeVisualTag {};
    struct C_CapsuleVisualTag {};
    struct C_CylinderVisualTag {};
    struct C_ConeVisualTag {};
    struct C_Collidable {};
    struct C_Radius { real_t r; };
    // Softbody visual surface: stores mesh topology (triangles) and the node entities driving vertices
    struct C_SoftBodyVisualTag {};
    struct C_SoftBodySurface {
        std::vector<Eigen::Vector3i> triangles;    // topology (indices into nodes)
        std::vector<entt::entity> nodes;           // one node per original OBJ vertex
    };
    // Beam element component
    struct C_BeamElement {
        std::optional<entt::entity> prev;
        std::optional<entt::entity> next;
        real_t l0{(real_t)0};
        real_t l{(real_t)0};
    };
    // Mesh components
    struct C_Mesh { std::string path; Vector3r scale{1,1,1}; };
    struct C_MeshVisualTag {};
    // HeightField components (static terrain)
    struct C_HeightField { std::string path; real_t x_dim{1}, y_dim{1}; real_t z_scale{1}; real_t min_height{0}; };
    struct C_HeightFieldVisualTag {};
    struct C_RB_HeightField { };
    // Rigid-body type components (exactly one per rigid body kind)
    struct C_RB_Cube { Vector3r center{Vector3r::Zero()}; Vector3r halfExtents; Quaternion4r q{Quaternion4r::Identity()}; };
    struct C_RB_Plane { Vector3r normal; Vector3r up; real_t sizeX; real_t sizeY; };
    struct C_RB_Mesh { };
    struct C_RB_Sphere { };
    struct C_RB_Capsule { real_t radius; real_t halfLength; };
    struct C_RB_Cylinder { real_t radius; real_t halfLength; };
    struct C_RB_Cone { real_t radius; real_t height; };

    // One-shot external wrench components (cleared every force rebuild)
    struct C_ExternalForce { Vector3r f; };
    struct C_ExternalTorque { Vector3r tau; };

    struct C_TrackTag {std::string name;};

    // Body index assigned by the assembler (stable across rebuilds unless structure changes)
    struct C_BodyIndex { int b; };

    // No numQ/numV or DOF accessors here; DynamicsAssembler owns DOF scans
    // Helper: number of dynamic bodies (with a body index)
    int numBodies() const;

    // State changed: positions/velocities modified
    void markStateDirty() const { m_state_dirty = true; }
    // Mark that structure changed (objects added/removed or dynamics tags changed)
    void markStructureDirty() const { m_structure_dirty = true; m_num_bodies_dirty = true; }
    // Forces changed: external forces like gravity updated
    void markForcesDirty() const { m_forces_dirty = true; }

    // Queries (non-consuming)
    bool isStateDirty() const { return m_state_dirty; }
    bool isStructureDirty() const { return m_structure_dirty; }
    bool isForcesDirty() const { return m_forces_dirty; }

    // Consumers (clear-on-read)
    bool consumeStateDirty() const { bool b = m_state_dirty; m_state_dirty = false; return b; }
    bool consumeStructureDirty() const { bool b = m_structure_dirty; m_structure_dirty = false; return b; }
    bool consumeForcesDirty() const { bool b = m_forces_dirty; m_forces_dirty = false; return b; }

    void applyForce(entt::entity e, const Vector3r& force_world, const Vector3r& torque_world);
    // Apply a pure moment specified in world coordinates; internally converted to body frame
    // (used when a torque is defined about an inertial-axis hinge).
    void applyInertialTorque(entt::entity e, const Vector3r& torque_world);
    void makeStatic(entt::entity e);

    // Minimal setters
    void setPosition(entt::entity e, const Vector3r& p);
    void setOrientation(entt::entity e, const Quaternion4r& q);
    void setLinearVelocity(entt::entity e, const Vector3r& v);
    void setAngularVelocity(entt::entity e, const Vector3r& w);
    void setVelocityByForce(entt::entity e, const Vector3r& v, const Vector3r& w);

    void explicitPositionUpdate(real_t dt);
    void linearImplicitPositionUpdate(real_t dt);

    void track(entt::entity e, const std::string& name);
    void writeTrackedStateToCsv(real_t t);

private:

    entt::registry m_reg;
    Vector3r m_gravity;  // gravity vector
    // no mass/structure caches here

    // System dirty flags (mutable so const getters/setters can flip)
    mutable bool m_state_dirty = true;     // q or v changed outside of physics loop
    mutable bool m_structure_dirty = true; // objects added/removed
    mutable bool m_forces_dirty = true;    // external forces changed

    // Cached number of bodies (entities with C_BodyIndex & C_PhysicsObject)
    mutable int m_num_bodies_cached = -1;
    mutable bool m_num_bodies_dirty = true;

    // assignDofs_ moved to DynamicsAssembler

    // Persistent subsystems
    config::Config m_cfg{}; // global config
    std::unique_ptr<collision::CollisionCoal> m_collision_mgr; // created on first use
    std::unique_ptr<cardillo::misc::TimingManager> m_timings;  // created on first use

    // Warmstart provider (strategy owned by system). Default implementation is WarmstartCache.
    std::unique_ptr<cardillo::solver::WarmstartProvider> m_warmstart_provider;
    // New constraint-pattern storage (owned by the system)
    std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>> m_constraints_new;

    // Asset manager (new abstraction)
    std::shared_ptr<class PhysicsAssets> m_assets;

public:
    void setAssets(std::shared_ptr<class PhysicsAssets> assets) { m_assets = std::move(assets); }
    class PhysicsAssets& assets();
    const class PhysicsAssets& assets() const;
};

} // namespace cardillo
