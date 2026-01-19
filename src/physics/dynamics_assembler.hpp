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
    explicit DynamicsAssembler(PhysicsSystem& sys) : m_sys(sys) {}

    // Legacy one-shot assembly helpers removed; use cached getters instead.

    // Cached getters (rebuild on demand) as concatenated global vectors
    const VectorXr& qVec();   // stacked positions
    const VectorXr& vVec();   // stacked velocities
    const VectorXr& fVec();   // stacked external forces
    const VectorXr& fVecExternal();   // stacked external forces (gravity + applied forces) only
    const VectorXr& fVecGyroscopic(); // stacked gyroscopic forces/torques only

    // Columns correspond to stacked body velocity DOFs in body order, laid out contiguously per body
    const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W() const { return m_W_sparse; }
    // Prefix sums of column offsets per body in stacked velocity vector; size = Nb + 1
    // Body b occupies columns [bodyVelOffsets()[b], bodyVelOffsets()[b+1])
    const std::vector<int>& bodyVelOffsets() const { return m_body_vel_offsets; }
    // Prefix sums for stacked position vector
    const std::vector<int>& bodyPosOffsets() const { return m_body_pos_offsets; }

    // Diagonal of block-diagonal Minv (size = total velocity dofs)
    const VectorXr& MinvDiag() const { return m_Minv_diag; }
    const VectorXr& MDiag() const { return m_M_diag; }
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
    void writePositionToSystem(const VectorXr& q);
    void writeVelocityToSystem(const VectorXr& v);
    void writeStateToSystem(const VectorXr& q, const VectorXr& v);
    // Assign contiguous DOF indices (3 per dynamic entity)
    void assignDofs();

    // Detect and set contacts from the current system state
    void updateContactsFromSystem();

    // Check system dirty flags and structural updates once per step and rebuild caches as needed
    void refreshState();
    void refreshCollisionsAndSprings(real_t dt, real_t theta, bool includeGyroInMatrix = false, bool lambdaTheta = false);
    void refreshCollisionsAndSpringsStormerVerlet(real_t dt);


private:
    PhysicsSystem& m_sys;

    // Cached data
    // Concatenated state
    VectorXr m_q_vec;
    VectorXr m_v_vec;
    VectorXr m_f_vec;
    VectorXr m_f_vec_external;   // gravity + external forces/torques only
    VectorXr m_f_vec_gyroscopic; // gyroscopic forces/torques only
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

   
    // Per-spring diagonal C = K^{-1} (scalar) (size = numSprings)
    VectorXr m_Cdiag; // size = numSprings
    // Per-damper diagonal A = D^{-1} (scalar) (size = numDampers)
    VectorXr m_Adiag; // size = numDampers

    // Extended block matrix S (sparse) and its sparse factorization
    CscMatrix m_S_sparse; // size extV x extV
    // Use SparseLU on the symmetric S matrix.
    std::optional<Eigen::SparseLU<CscMatrix>> m_S_sparse_lu; // lazily constructed

    // Store Lagrange multipliers (they are being integrated)
    VectorXr m_Lambda_g;
    VectorXr m_Lambda_gamma;

public:
    // Accessors for newly added matrices
    const Eigen::SparseMatrix<real_t>& WgSparse() const { return m_Wg; }
    const Eigen::SparseMatrix<real_t>& WgammaSparse() const { return m_Wgamma; }

    // Per-spring diagonal C/A
    const VectorXr& Cdiag() const { return m_Cdiag; }
    const VectorXr& Adiag() const { return m_Adiag; }

    // Counts
    // Number of spring rows (rows in m_Wg / length of m_Cdiag)
    int numSprings() const { return (int)m_Cdiag.size(); }
    int numDampers() const { return (int)m_Adiag.size(); }

    // Accessors for Lagrange multipliers
    const VectorXr& Lambda_g() const { return m_Lambda_g; }
    const VectorXr& Lambda_gamma() const { return m_Lambda_gamma; }
    void setLambda_g(const VectorXr& lam);
    void setLambda_gamma(const VectorXr& lam) { m_Lambda_gamma = lam; }

    // Extended block matrix S (sparse) accessor and solver
    const CscMatrix& SSparse() const { return m_S_sparse; }

    // Build and factorize the effective mass matrix S = M + dt^2 * Wg * K * Wg^T + h * W_gamma * D * W_gamma^T
    // Returns true on successful factorization.
    bool buildAndFactorS(real_t dt, real_t theta, bool includeGyroInMatrix = false, bool lambdaTheta = false);
    // Solve the full extended system S * x = rhs_ext and return the complete solution (ext-length)
    VectorXr solveS(const VectorXr& rhs_ext) const;

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
    void rebuildInteractionW_();
    bool buildAndFactorS_StormerVerlet(real_t dt);
};

} // namespace cardillo::physics
