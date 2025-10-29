#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <optional>
#include <functional>
#include <unordered_map>
#include <entt/entt.hpp>
#include <petscsys.h>
#include "../misc/types.hpp"
#include "../misc/dofs.hpp"
#include "../config/config.hpp"
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

namespace cardillo {

// A minimal, standard-C++ physics system for frictionless point masses
// with translational DOFs only (no rotations).
class PhysicsSystem {
public:
    PhysicsSystem();
    explicit PhysicsSystem(const config::Config& cfg);
    ~PhysicsSystem();
    // Global config accessible across subsystems
    void setConfig(const config::Config& cfg) { m_cfg = cfg; setGravity(m_cfg.sim_gravity); }
    const config::Config& config() const { return m_cfg; }

    // Persistent collision manager (COAL) storage and access
    collision::CollisionCoal& collisionManager();
    const collision::CollisionCoal& collisionManager() const;


    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const { return m_gravity; }

    // Visual/collision plane configuration (for future use in collisions)
    struct Plane {
        Vector3r center{0,0,0};
        Vector3r normal{0,0,1};
        Vector3r up{0,1,0};
        real_t sizeX{5}, sizeY{5}; // half extents for visualization
    };
    // Axis-aligned cube visual with quaternion orientation
    struct Cube {
        Vector3r center{0,0,0};
        Vector3r halfExtents{0.5,0.5,0.5};
    Quaternion4r q = Quaternion4r::Identity(); // orientation
    };
    struct Capsule {
        real_t radius{0.5};
        real_t halfLength{0.5};
    };
    // Obstacle (static) visuals
    index_t addObstacleBody(const Plane& p);
    index_t addObstacleBody(const Cube& c);
    // Static mesh obstacle (visual + collider), no dynamics/inertia normalization
    index_t addObstacleMesh(const Vector3r& position,
                            const Quaternion4r& orientation,
                            const std::string& meshPath,
                            const Vector3r& scale = Vector3r(1,1,1));
    // Dynamic rigid body creation: pose and spatial velocity
    index_t addRigidBody(real_t mass,
                         const Vector3r& position,
                         const Quaternion4r& orientation,
                         const Vector3r& linearVelocity,
                         const Vector3r& angularVelocity,
                         const Cube& shape);
    index_t addRigidBodyCapsule(real_t mass,
                                const Vector3r& position,
                                const Quaternion4r& orientation,
                                const Vector3r& linearVelocity,
                                const Vector3r& angularVelocity,
                                const Capsule& shape);
    // Mesh-based rigid body (collision via COAL BVH, visual via VTK mesh output)
    index_t addRigidBodyMesh(real_t mass,
                             const Vector3r& position,
                             const Quaternion4r& orientation,
                             const Vector3r& linearVelocity,
                             const Vector3r& angularVelocity,
                             const std::string& meshPath,
                             const Vector3r& scale = Vector3r(1,1,1));
    // Static HeightField obstacle sourced from an EXR heightmap
    index_t addObstacleHeightField(const Vector3r& position,
                                   const Quaternion4r& orientation,
                                   const std::string& exrPath,
                                   real_t x_dim,
                                   real_t y_dim,
                                   real_t z_scale = (real_t)1.0,
                                   real_t min_height = (real_t)0.0);
    // Sphere-based rigid body (with rotation and inertia)
    index_t addRigidBodySphere(real_t mass,
                               const Vector3r& position,
                               const Quaternion4r& orientation,
                               const Vector3r& linearVelocity,
                               const Vector3r& angularVelocity,
                               real_t radius);
  
    index_t addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius = (real_t)0.05);

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
    VectorXr getForce( entt::entity e ) const;        // Linear and angular combined

    // Shared mesh asset access (deduplicates loading across VTK/collision/inertia)
    struct MeshAsset {
        coal::BVHModelPtr_t bvh;                    // scaled geometry BVH
        std::vector<Eigen::Vector2f> uvs;           // optional per-vertex UVs
        bool hasUV = false;
        // Normalization metadata (unit density)
        Vector3r inertia_diag_unit = Vector3r::Zero(); // principal inertias for rho=1, about COM
        real_t volume = (real_t)0.0;                    // signed volume (>0 for outward-facing)
        // Additional normalization info
        Matrix33r Rpa = Matrix33r::Identity();          // principal-axes rotation used when normalized
        Vector3r com = Vector3r::Zero();                // center of mass in scaled mesh frame
        bool normalized = false;                        // whether mesh vertices in bvh are normalized (COM-centered, PA-aligned)
    };
    // Entity-based accessor: determines static/dynamic and may adjust pose for dynamic meshes
    const MeshAsset& getMeshAsset(entt::entity e) const;
    // Shared heightfield asset (deduplicates loading from EXR)
    struct HeightFieldAsset {
        std::shared_ptr<coal::HeightField<coal::AABB>> hf; // collider geometry
        int rows{0}, cols{0};
        real_t x_dim{1}, y_dim{1};
        real_t z_scale{1}, min_height{0};
        std::string path;
    };
    const HeightFieldAsset& getHeightFieldAsset(entt::entity e) const;

    // Optional inertia component for rigid bodies with arbitrary shapes
    struct C_InertiaDiag { Vector3r I; }; // body-frame diagonal inertia (Ixx,Iyy,Izz)

    // Expose ECS for external querying (read-only)
    const entt::registry& ecs() const { return m_reg; }

    // Access to warmstart provider owned by the system (may be nullptr)
    cardillo::solver::WarmstartProvider* warmstartProvider() const { return m_warmstart_provider.get(); }

    // Public ECS component/tag types for queries
    struct C_Mass { real_t m; };
    struct C_Position3 { Vector3r value; };
    struct C_LinearVelocity3 { Vector3r value; };
    struct C_AngularVelocity3 { Vector3r value; };
    struct C_Orientation { Quaternion4r q; };
    struct C_PhysicsObject {};
    struct C_PointMassTag {};
    struct C_RigidBodyTag {};
    struct C_Plane { Vector3r normal; Vector3r up; real_t sizeX; real_t sizeY; };
    struct C_Cube { Vector3r halfExtents; };
    struct C_Capsule { real_t radius; real_t halfLength; };
    struct C_Friction { real_t mu; }; // optional friction coefficient per entity (>=0), absent => 0
    struct C_VisualObject {};
    struct C_PointVisualTag {};
    struct C_PlaneVisualTag {};
    struct C_CubeVisualTag {};
    struct C_CapsuleVisualTag {};
    struct C_Collidable {};
    struct C_Radius { real_t r; };
    // Mesh components
    struct C_Mesh { std::string path; Vector3r scale{1,1,1}; };
    struct C_MeshVisualTag {};
    // HeightField components (static terrain)
    struct C_HeightField { std::string path; real_t x_dim{1}, y_dim{1}; real_t z_scale{1}; real_t min_height{0}; };
    struct C_HeightFieldVisualTag {};
    struct C_RB_HeightField { };
    // Rigid-body type components (exactly one per rigid body kind)
    struct C_RB_Cube { Vector3r halfExtents; };
    struct C_RB_Plane { Vector3r normal; Vector3r up; real_t sizeX; real_t sizeY; };
    struct C_RB_Mesh { };
    struct C_RB_Sphere { };
    struct C_RB_Capsule { real_t radius; real_t halfLength; };

    // Body index assigned by the assembler (stable across rebuilds unless structure changes)
    struct C_BodyIndex { int b; };

    // No numQ/numV or DOF accessors here; DynamicsAssembler owns DOF scans
    // Helper: number of dynamic bodies (with a body index)
    int numBodies() const;

    // State changed: positions/velocities modified
    void markStateDirty() const { m_state_dirty = true; }
    // Structure changes are picked up by assemblers directly from ECS; no flag needed
    void markStructureDirty() const { /* no-op, kept for API compatibility */ }
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

private:
    // Helper to add common rigid-body components
    void emplaceRigidBodyCommon_(entt::entity e,
                                 real_t mass,
                                 const Vector3r& position,
                                 const Quaternion4r& orientation,
                                 const Vector3r& linearVelocity,
                                 const Vector3r& angularVelocity);

private:

    entt::registry m_reg;
    Vector3r m_gravity;  // gravity vector
    // no mass/structure caches here

    // System dirty flags (mutable so const getters/setters can flip)
    mutable bool m_state_dirty = true;     // q or v changed outside of physics loop
    mutable bool m_structure_dirty = true; // objects added/removed
    mutable bool m_forces_dirty = true;    // external forces changed

    // assignDofs_ moved to DynamicsAssembler
    entt::entity createRigidVisualEntity_(const Vector3r& center);

    // Persistent subsystems
    config::Config m_cfg{}; // global config
    std::unique_ptr<collision::CollisionCoal> m_collision_mgr; // created on first use

    // Warmstart provider (strategy owned by system). Default implementation is WarmstartCache.
    std::unique_ptr<cardillo::solver::WarmstartProvider> m_warmstart_provider;

    // Shared mesh cache (keyed by path + scale)
    mutable std::unordered_map<std::string, MeshAsset> m_meshCache;
    // Shared heightfield cache (keyed by path + dims + scales)
    mutable std::unordered_map<std::string, HeightFieldAsset> m_hfCache;
};

} // namespace cardillo
