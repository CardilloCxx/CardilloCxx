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
#include "misc/types.hpp"
#include "api/physics_types.hpp"
#include "ecs_types.hpp"
#include "../misc/spline.hpp"
#include "assets/assets.hpp"
#include "../misc/dofs.hpp"
#include "../config/config.hpp"
#include "../misc/timings/TimingManager.hpp"
#include <coal/BVH/BVH_model.h>
#include <coal/collision_object.h>
#include <coal/broadphase/broadphase.h>
#include <coal/hfield.h>

namespace cardillo { namespace solver { class WarmstartProvider; } }
namespace cardillo { namespace collision { class CollisionCoal; } }
namespace cardillo { namespace physics { class ConstraintPattern; class LinearDistanceConstraint; struct JointFrame; } }

namespace cardillo {
class World {
public:
    static Quaternion4r alignQuaternionTo(const Quaternion4r& q_in, const Quaternion4r& q_ref) {
        Quaternion4r q = q_in;
        if (!q.coeffs().allFinite()) return q_ref;
        q.normalize();
        if (q_ref.dot(q) < (real_t)0) q.coeffs() = -q.coeffs();
        return q;
    }

    using RigidState = cardillo::physics::RigidState;
    using RigidShape = cardillo::physics::RigidShape;
    using RigidProps = cardillo::physics::RigidProps;
    using BeamCrossSection = cardillo::physics::BeamCrossSection;
    using BeamSpringParams = cardillo::physics::BeamSpringParams;
    using BeamBodyType = cardillo::physics::BeamBodyType;

    World();
    explicit World(const config::Config& cfg);
    ~World();
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

    // Non-owning external pointers (set when an external manager owns the subsystems)
    collision::CollisionCoal* m_collision_mgr_external{nullptr};
    cardillo::misc::TimingManager* m_timings_external{nullptr};
    cardillo::solver::WarmstartProvider* m_warmstart_provider_external{nullptr};

public:
    void setAssets(std::shared_ptr<class PhysicsAssets> assets) { m_assets = std::move(assets); }
    class PhysicsAssets& assets();
    const class PhysicsAssets& assets() const;

    // Set non-owning pointers to subsystems owned by an external lifecycle manager (e.g., PhysicsEngine)
    void setCollisionManager(collision::CollisionCoal* mgr);
    void setTimings(cardillo::misc::TimingManager* timings);
    void setWarmstartProvider(cardillo::solver::WarmstartProvider* provider);
};

} // namespace cardillo
