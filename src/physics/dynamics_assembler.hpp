#pragma once

#include <optional>
#include <vector>
#include <Eigen/SparseCore>

#include "../misc/types.hpp"
#include "physics_system.hpp"
#include "../collision/types.hpp"

namespace cardillo::physics {


// Contracts:
// - Provides cached assembly of state vectors (q, v), force vector f, mass matrix M,
//   contact matrix W (maps v -> normal velocities), and Delassus matrix G = W M^{-1} W^T.
// - Maintains dirty flags; callers can mark state/structure/forces/contacts dirty or set contacts.
// - Rebuilds lazily on demand.
class DynamicsAssembler {
public:
    explicit DynamicsAssembler(const PhysicsSystem& sys) : m_sys(sys) {}

    // Legacy one-shot assembly helpers removed; use cached getters instead.

    // Cached getters (rebuild on demand) as concatenated global vectors
    const VectorXr& qVec();   // stacked positions
    const VectorXr& vVec();   // stacked velocities
    const VectorXr& fVec();   // stacked external forces

    // Columns correspond to stacked body velocity DOFs in body order, laid out contiguously per body
    const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W() const { return m_W_sparse; }
    // Prefix sums of column offsets per body in stacked velocity vector; size = Nb + 1
    // Body b occupies columns [bodyVelOffsets()[b], bodyVelOffsets()[b+1])
    const std::vector<int>& bodyVelOffsets() const { return m_body_vel_offsets; }
    // Prefix sums for stacked position vector
    const std::vector<int>& bodyPosOffsets() const { return m_body_pos_offsets; }

    // Diagonal of block-diagonal Minv (size = total velocity dofs)
    const VectorXr& MinvDiag() const { return m_Minv_diag; }
    // Map dynamic-only contact index (used in W/G) back to original contact index in m_contacts
    const std::vector<int>& dynamicContactToOriginalAll() const { return m_contact_index_orig; }
    // Access underlying system (for debug / diagnostics)
    const PhysicsSystem& system() const { return m_sys; }
    // Expose current contacts (includes penetration and points) for biasing, debug, etc.
    const std::vector<collision::Contact>& contacts() const { return m_contacts; }

    // DOF queries owned by the assembler (scan ECS indices)
    index_t numQ() const { return m_numQ; }
    index_t numV() const { return m_numV; }
    
    // Explicit state IO
    void loadStateFromSystem();
    void writeStateToSystem(const VectorXr& q_concat, const VectorXr& v_concat);
    // Assign contiguous DOF indices (3 per dynamic entity)
    void assignDofs();

    // Detect and set contacts from the current system state
    void updateContactsFromSystem();

    // Check system dirty flags and structural updates once per step and rebuild caches as needed
    void refreshState();
    void refreshCollisionsAndSprings( real_t dt );


private:
    const PhysicsSystem& m_sys;

    // Cached data
    // Concatenated state
    VectorXr m_q_vec;
    VectorXr m_v_vec;
    VectorXr m_f_vec;
    // Legacy caches removed
    std::vector<VectorXr> m_v_compat; // kept temporarily for transition; unused in solvers

    // Sparse contact Jacobian and supporting mappings
    Eigen::SparseMatrix<real_t, Eigen::RowMajor> m_W_sparse; // (C_dynamic x totalV), RowMajor to allow efficient row access
    std::vector<int> m_body_vel_offsets;    // size Nb+1, prefix sums for body velocity columns
    std::vector<int> m_body_pos_offsets;    // size Nb+1, prefix sums for body position columns
    VectorXr m_Minv_diag;                   // size totalV, diagonal of Minv
    VectorXr m_M_diag;                      // size totalV, diagonal of M (non-inverted)

    // Additional per-contact/block matrices used by solvers
    // Wg maps constraint-space generalized g/constraint DOFs to global velocity DOFs (same shape as m_W_sparse)
    Eigen::SparseMatrix<real_t> m_Wg;
    // W_gamma often is the gamma-related Jacobian (negated Wg for one side); keep as sparse
    Eigen::SparseMatrix<real_t> m_Wgamma;

    // Contact-space stiffness/constraint matrices (sparse). Dimensions: C_dyn x C_dyn
    Eigen::SparseMatrix<real_t> m_K; // stiffness in constraint space (C_dyn x C_dyn)
    Eigen::SparseMatrix<real_t> m_D; // damping in constraint space (C_dyn x C_dyn)

    // Effective mass S in velocity DOF space and its Cholesky factorization
    Eigen::SparseMatrix<real_t> m_S; // size totalV x totalV
    std::optional<Eigen::SimplicialLLT<Eigen::SparseMatrix<real_t>>> m_S_llt; // lazily constructed

    // Concatenated spring deformations and rates
    VectorXr m_gcat;
    VectorXr m_gdotcat;

public:
    // Accessors for newly added matrices
    const Eigen::SparseMatrix<real_t>& WgSparse() const { return m_Wg; }
    const Eigen::SparseMatrix<real_t>& WgammaSparse() const { return m_Wgamma; }
    const Eigen::SparseMatrix<real_t>& Kmat() const { return m_K; }
    const Eigen::SparseMatrix<real_t>& Dmat() const { return m_D; }

    // Accessors for concatenated spring deformations and rates
    const VectorXr& gcat() const { return m_gcat; }
    const VectorXr& gdotcat() const { return m_gdotcat; }

    // Build and factorize the effective mass matrix S = M + dt^2 * Wg * K * Wg^T + h * W_gamma * D * W_gamma^T
    // Returns true on successful factorization.
    bool buildAndFactorS(real_t dt);
    // Solve S x = rhs using the precomputed factorization. Factorization must be valid.
    VectorXr solveWithS(const VectorXr& rhs) const;

    std::vector<int> m_contact_index_orig; // size Nc_dynamic, maps dynamic contact id -> original contact index
    std::vector<collision::Contact> m_contacts;

    // Dirty flags
    // (unused) bool m_contacts_dirty{true};

    // Cached sizes
    index_t m_numQ{0};
    index_t m_numV{0};

    // Rebuild helpers
    void rebuildMass_();
    void rebuildForces_();
    void rebuildW_();
    void rebuidInteractionW_();
};

} // namespace cardillo::physics
