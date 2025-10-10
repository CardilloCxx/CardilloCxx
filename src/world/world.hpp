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
// - midpoint-rule time integration under uniform gravity
class PhysicsSystem {
public:
    struct PointMass {
        real_t m;           // mass
        Vector3r x;         // position
        Vector3r v;         // velocity
        // DoF indices into global vectors (3 translational per mass)
        index_t q_start = -1;
        index_t v_start = -1; // identical to q_start in this translational-only model
    };

    // Create an empty system
    PhysicsSystem();

    // Set gravity (default: [0,0,-9.81])
    void setGravity(const Vector3r& g);

    // Add a point mass; returns its index
    index_t addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0);

    // Number of translational DoFs (= 3 * numPointMasses)
    index_t numQ() const { return static_cast<index_t>(m_q_dofs); }
    index_t numV() const { return static_cast<index_t>(m_v_dofs); }

    // Assemble diagonal mass matrix M (size: VxV)
    Eigen::SparseMatrix<real_t> assembleMassMatrix() const;

    // Assemble generalized force vector f (size: V)
    // Currently only gravity: f_i = m_i * g
    VectorXr assembleForceVector() const;

    // Getters for state packing/unpacking
    VectorXr packQ() const; // positions stacked [x1;y1;z1; x2;...]
    VectorXr packV() const; // velocities stacked
    void unpackQ(const RefVectorXr& q);
    void unpackV(const RefVectorXr& v);

    // Midpoint rule integration for dt
    // q_{n+1} = q_n + dt * v_mid
    // v_{n+1} = v_n + dt * a_mid,  with a_mid = M^{-1} f(q_mid, v_mid)
    void stepMidpoint(real_t dt);

    // Access to individual masses if needed
    const std::vector<PointMass>& masses() const { return m_masses; }

private:
    std::vector<PointMass> m_masses;
    Vector3r m_gravity;  // gravity vector
    size_t m_q_dofs = 0;
    size_t m_v_dofs = 0;

    // Internal helpers
    void assignDofs_();
};

} // namespace cardillo