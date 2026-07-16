#pragma once

#include <coal/BVH/BVH_model.h>
#include <coal/broadphase/broadphase.h>
#include <coal/collision_object.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <entt/entt.hpp>
#include <functional>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>
#include "../config/config.hpp"
#include "../misc/dofs.hpp"
#include "../misc/math_helper.hpp"
#include "../misc/spline.hpp"
#include "../misc/timings/TimingManager.hpp"
#include "api/physics_types.hpp"
#include "assets/assets.hpp"
#include "ecs_types.hpp"
#include "misc/types.hpp"

namespace cardillo {
namespace solver {
class WarmstartProvider;
}
}  // namespace cardillo
namespace cardillo {
namespace collision {
class CollisionCoal;
}
}  // namespace cardillo
namespace cardillo {
namespace physics {
class ConstraintPattern;
class LinearDistanceConstraint;
struct JointFrame;
}  // namespace physics
}  // namespace cardillo

namespace cardillo {
class World {
   public:
    explicit World(const config::Config& cfg);
    ~World();

    entt::registry m_reg;
    Vector3r m_gravity;  // gravity vector

    // System dirty flags (mutable so const getters/setters can flip)
    mutable bool m_state_dirty = true;      // q or v changed outside of physics loop
    mutable bool m_structure_dirty = true;  // objects added/removed
    mutable bool m_forces_dirty = true;     // external forces changed

    // Cached number of bodies (entities with C_BodyIndex & C_PhysicsObject)
    mutable int m_num_bodies_cached = -1;
    mutable bool m_num_bodies_dirty = true;

    config::Config m_cfg{};  // global config

    std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>> m_constraints_new;
    std::shared_ptr<class PhysicsAssets> m_assets;

    // Access global config
    const config::Config& config() const { return m_cfg; }
    config::Config& config() { return m_cfg; }

    // Access asset manager
    class PhysicsAssets& assets();
    const class PhysicsAssets& assets() const;

    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const { return m_gravity; }

    // Dynamics getters (Cache them inside the entity to avoid recomputation)
    // Mass diagonal (vector) for the entity's velocity-space dofs.
    // Order matches getVelocity():
    //  - Rigid body: [m, m, m, Ixx, Iyy, Izz]
    //  - Point mass: [m, m, m]
    VectorXr getMassDiag(entt::entity e) const;
    MatrixXXr getMass(entt::entity e) const;  // Linear Inertia and Angular Inertia
    // Inverse mass diagonal (vector), see getMassDiag
    VectorXr getMassInverseDiag(entt::entity e) const;
    // Generalized inertia diagonal getter (body frame). Computes or retrieves best-known diag.
    Vector3r getInertiaDiag(entt::entity e) const;
    VectorXr getPosition(entt::entity e) const;         // Linear and angular combined
    VectorXr getVelocity(entt::entity e) const;         // Linear and angular combined
    VectorXr getForce(entt::entity e) const;            // Linear and angular combined (gravity + external + gyroscopic)
    VectorXr getForceExternal(entt::entity e) const;    // Gravity + external forces/torques only
    VectorXr getForceGyroscopic(entt::entity e) const;  // Gyroscopic torque only (tau = -w x (I*w))
    real_t getKineticEnergy(entt::entity e) const;

    // Shared asset access (wrappers over PhysicsAssets using entity components)
    const ::cardillo::MeshAsset& getMeshAsset(entt::entity e) const;

    // Expose ECS for external querying (read-only)
    const entt::registry& ecs() const { return m_reg; }
    // Mutable ECS access when external systems need to emplace components
    entt::registry& ecs() { return m_reg; }

    // Access all constraint patterns (mutable and const)
    std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>>& constraintPatterns() { return m_constraints_new; }
    const std::vector<std::unique_ptr<cardillo::physics::ConstraintPattern>>& constraintPatterns() const { return m_constraints_new; }

    int numBodies() const;

    // State changed: positions/velocities modified
    void markStateDirty() const { m_state_dirty = true; }
    // Mark that structure changed (objects added/removed or dynamics tags changed)
    void markStructureDirty() const {
        m_structure_dirty = true;
        m_num_bodies_dirty = true;
    }
    // Forces changed: external forces like gravity updated
    void markForcesDirty() const { m_forces_dirty = true; }

    // Queries (non-consuming)
    bool isStateDirty() const { return m_state_dirty; }
    bool isStructureDirty() const { return m_structure_dirty; }
    bool isForcesDirty() const { return m_forces_dirty; }

    // Consumers (clear-on-read)
    bool consumeStateDirty() const {
        bool b = m_state_dirty;
        m_state_dirty = false;
        return b;
    }
    bool consumeStructureDirty() const {
        bool b = m_structure_dirty;
        m_structure_dirty = false;
        return b;
    }
    bool consumeForcesDirty() const {
        bool b = m_forces_dirty;
        m_forces_dirty = false;
        return b;
    }

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
    void setTrajectory(entt::entity e, std::optional<std::function<TrajectoryPose(real_t)>> positionFunc, std::optional<std::function<TrajectoryTwist(real_t)>> velocityFunc);
    void removeTrajectory(entt::entity e);

    void track(entt::entity e, const std::string& name);
};

}  // namespace cardillo
