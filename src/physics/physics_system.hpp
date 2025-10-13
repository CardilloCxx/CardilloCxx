#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <optional>
#include <functional>
#include <entt/entt.hpp>
#include <petscsys.h>
#include "../misc/types.hpp"
#include "../misc/dofs.hpp"

namespace cardillo {

// A minimal, standard-C++ physics system for frictionless point masses
// with translational DOFs only (no rotations).
class PhysicsSystem {
public:
    PhysicsSystem();

    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const { return m_gravity; }

    // Visual/collision plane configuration (for future use in collisions)
    struct Plane {
        Vector3r center{0,0,0};
        Vector3r normal{0,0,1};
        Vector3r up{0,1,0};
        real_t sizeX{5}, sizeY{5}; // half extents for visualization
    };
    // Axis-aligned cube visual
    struct Cube {
        Vector3r center{0,0,0};
        Vector3r halfExtents{0.5,0.5,0.5};
        Matrix33r R = Matrix33r::Identity(); // orientation (columns are local axes)
    };
    // Rigid visuals (no physics yet)
    index_t addRigidBody(const Plane& p);
    index_t addRigidBody(const Cube& c);
    index_t addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0, real_t radius = (real_t)0.05);

    // Expose ECS for external querying (read-only)
    const entt::registry& ecs() const { return m_reg; }

    // Public ECS component/tag types for queries
    struct C_Mass { real_t m; };
    struct C_Position3 { Vector3r value; };
    struct C_LinearVelocity3 { Vector3r value; };
    struct C_PhysicsObject {};
    struct C_PointMassTag {};
    struct C_RigidBodyTag {};
    struct C_Plane { Vector3r normal; Vector3r up; real_t sizeX; real_t sizeY; };
    struct C_Cube { Vector3r halfExtents; Matrix33r R; };
    struct C_VisualObject {};
    struct C_PointVisualTag {};
    struct C_PlaneVisualTag {};
    struct C_CubeVisualTag {};
    struct C_Collidable {};
    struct C_Radius { real_t r; };
    // DOF indices as components (public so assemblers can read them)
    struct C_PositionIndex3 { Index<3> idx; };
    struct C_LinearVelocityIndex3 { Index<3> idx; };

    // No numQ/numV or DOF accessors here; DynamicsAssembler owns DOF scans

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

    entt::registry m_reg;
    Vector3r m_gravity;  // gravity vector
    // no mass/structure caches here

    // System dirty flags (mutable so const getters/setters can flip)
    mutable bool m_state_dirty = true;     // q or v changed outside of physics loop
    mutable bool m_structure_dirty = true; // objects added/removed
    mutable bool m_forces_dirty = true;    // external forces changed

    // assignDofs_ moved to DynamicsAssembler
    entt::entity createRigidVisualEntity_(const Vector3r& center);
};

} // namespace cardillo
