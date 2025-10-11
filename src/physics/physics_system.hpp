#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <functional>
#include <entt/entt.hpp>
#include <petscsys.h>
#include "../misc/types.hpp"
#include "../misc/dofs.hpp"

namespace cardillo {

// A minimal, standard-C++ physics system for frictionless point masses
// with translational DOFs only (no rotations). Supports:
// - adding point masses
// - assembling the (diagonal) mass matrix
// - providing pack/unpack helpers for solvers
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
    };
    // Rigid visuals (no physics yet)
    index_t addRigidBody(const Plane& p);
    index_t addRigidBody(const Cube& c);
    void clearPlanes() { /* legacy no-op */ }
    bool hasPlanes() const { return false; }
    const std::vector<Plane> planes() const { return {}; }

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
    struct C_Cube { Vector3r halfExtents; };
    struct C_VisualObject {};
    struct C_PointVisualTag {};
    struct C_PlaneVisualTag {};
    struct C_CubeVisualTag {};

    // Add a point mass; returns an entity id encoded as index_t (for now)
    index_t addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0);

    index_t numQ() const { return m_q_dofs; }
    index_t numV() const { return m_v_dofs; }

    // Assemble diagonal mass matrix M (size: VxV)
    Eigen::SparseMatrix<real_t> assembleMassMatrix() const;

    // Assemble generalized force vector f (size: V); currently only gravity
    VectorXr assembleForceVector() const;

    // Cached mass diagonal (3*m per mass). Recomputed when masses/DOFs change.
    const VectorXr& massDiagonal() const;

    // State helpers
    VectorXr packQ() const; // positions stacked
    VectorXr packV() const; // velocities stacked
    void unpackQ(const RefVectorXr& q);
    void unpackV(const RefVectorXr& v);

    // Removed public snapshot masses(); writer should query via forEachVisualPoint.

private:
    // DOF indices as components (internal)
    // DOF indices as components
    struct C_PositionIndex3 { Index<3> idx; };
    struct C_LinearVelocityIndex3 { Index<3> idx; };

    entt::registry m_reg;
    Vector3r m_gravity;  // gravity vector
    index_t m_q_dofs = 0;
    index_t m_v_dofs = 0;
    mutable VectorXr m_Mdiag;     // cached mass diagonal
    mutable bool m_mass_dirty = true;

    void assignDofs_();
    entt::entity createRigidVisualEntity_(const Vector3r& center);
};

} // namespace cardillo
