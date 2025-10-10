#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "../misc/types.hpp"

namespace cardillo {

// A minimal, standard-C++ physics system for frictionless point masses
// with translational DOFs only (no rotations). Supports:
// - adding point masses
// - assembling the (diagonal) mass matrix
// - providing pack/unpack helpers for solvers
class PhysicsSystem {
public:
    struct PointMass {
        real_t m;           // mass
        Vector3r x;         // position
        Vector3r v;         // velocity
        index_t q_start = -1; // position DOF start
        index_t v_start = -1; // velocity DOF start (same as q_start here)
    };

    PhysicsSystem();

    void setGravity(const Vector3r& g);
    const Vector3r& gravity() const { return m_gravity; }

    // Add a point mass; returns its index
    index_t addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0);

    index_t numQ() const { return static_cast<index_t>(m_q_dofs); }
    index_t numV() const { return static_cast<index_t>(m_v_dofs); }

    // Assemble diagonal mass matrix M (size: VxV)
    Eigen::SparseMatrix<real_t> assembleMassMatrix() const;

    // Assemble generalized force vector f (size: V); currently only gravity
    VectorXr assembleForceVector() const;

    // State helpers
    VectorXr packQ() const; // positions stacked
    VectorXr packV() const; // velocities stacked
    void unpackQ(const RefVectorXr& q);
    void unpackV(const RefVectorXr& v);

    // Access to individual masses
    const std::vector<PointMass>& masses() const { return m_masses; }
    std::vector<PointMass>& masses() { return m_masses; }

private:
    std::vector<PointMass> m_masses;
    Vector3r m_gravity;  // gravity vector
    size_t m_q_dofs = 0;
    size_t m_v_dofs = 0;

    void assignDofs_();
};

} // namespace cardillo
