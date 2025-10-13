#pragma once

#include <optional>
#include <vector>
#include <Eigen/SparseCore>

#include "../misc/types.hpp"
#include "physics_system.hpp"
#include "../collision/collision.hpp"

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

    // Cached getters (rebuild on demand)
    const VectorXr& q();
    const VectorXr& v();
    const VectorXr& f();
    const Eigen::SparseMatrix<real_t>& M();
    const Eigen::SparseMatrix<real_t>& Minv();
    const VectorXr& MinvDiag();
    const Eigen::SparseMatrix<real_t>& W();
    const Eigen::SparseMatrix<real_t>& G();

    // DOF queries owned by the assembler (scan ECS indices)
    index_t numQ() const { return m_numQ; }
    index_t numV() const { return m_numV; }
    std::optional<index_t> velDofStart(entt::entity e) const;

    // Update inputs / invalidate caches
    void setContacts(std::vector<collision::Contact> contacts);
    // Dirty markers: contacts (others live on PhysicsSystem)
    void markContactsDirty();   // contact set or normals changed

    // Explicit state IO
    void loadStateFromSystem();
    void writeStateToSystem(const RefVectorXr& q, const RefVectorXr& v);

    // Assign contiguous DOF indices (3 per dynamic entity)
    void assignDofs();

    // Detect and set contacts from the current system state
    void updateContactsFromSystem();

    // Check system dirty flags and structural updates once per step and rebuild caches as needed
    void refreshState();

private:
    const PhysicsSystem& m_sys;

    // Cached data
    VectorXr m_q;
    VectorXr m_v;
    VectorXr m_f;
    Eigen::SparseMatrix<real_t> m_M;
    Eigen::SparseMatrix<real_t> m_Minv; // diagonal inverse mass as sparse
    VectorXr m_MinvDiag;                 // inverse mass diagonal as vector
    VectorXr m_Mdiag;                    // mass diagonal as vector
    Eigen::SparseMatrix<real_t> m_W;    // size C x V
    Eigen::SparseMatrix<real_t> m_G;    // size C x C
    std::vector<collision::Contact> m_contacts;

    // Dirty flags
    bool m_contacts_dirty{true};
    bool m_G_dirty{true};

    // Cached sizes
    index_t m_numQ{0};
    index_t m_numV{0};

    // Rebuild helpers
    void rebuildMass_();
    void rebuildForces_();
    void rebuildW_();
    void rebuildG_();
};

} // namespace cardillo::physics
