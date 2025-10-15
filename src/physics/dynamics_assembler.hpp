#pragma once

#include <optional>
#include <vector>
#include <Eigen/SparseCore>

#include "../misc/types.hpp"
#include "physics_system.hpp"
#include "../collision/types.hpp"

namespace cardillo::physics {

struct WBlockEdge {
    int bodyA;
    int bodyB;
    MatrixXXr WblockA; // dofP x dofV (currently 1x3)
    MatrixXXr WblockB; // dofP x dofV (currently 1x3)
};

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
    const std::vector<VectorXr>& q();
    const std::vector<VectorXr>& v();
    const std::vector<VectorXr>& f();

    // Block-level accessors (per-entity velocity blocks)
    // Current implementation assumes 3 dofs per entity; structure can be generalized later.
    const std::vector<VectorXr>& vBlocks();          // size Nb, each VectorXr size dofV
    const std::vector<MatrixXXr>& MinvBlocks();      // size Nb, each MatrixXXr dofVxdofV
    const std::vector<MatrixXXr>& WBlocks();         // size Nc * 2 each MatrixXXr dofPxdofV
    const std::pair<int, int> & WBlocksFromContact(index_t contact);  
    const std::vector<int>& WBlocksFromBody (index_t body); 
    // Bulk accessors to internal mappings (avoid per-element calls in tight loops)
    const std::vector<std::pair<int,int>>& WBlocksFromContactAll() const { return m_W_from_contact; }
    const std::vector<std::vector<int>>& WBlocksFromBodyAll() const { return m_W_from_body; }
    const std::vector<int>& WBlockToBodyAll() const { return m_W_block_to_body; }
    const std::vector<int>& WBlockToContactAll() const { return m_W_block_to_contact; }

    // DOF queries owned by the assembler (scan ECS indices)
    index_t numQ() const { return m_numQ; }
    index_t numV() const { return m_numV; }

    // Update inputs / invalidate caches
    void setContacts(std::vector<collision::Contact> contacts);
    // Dirty markers: contacts (others live on PhysicsSystem)
    void markContactsDirty();   // contact set or normals changed

    // Explicit state IO
    void loadStateFromSystem();
    void writeStateToSystem(const std::vector<VectorXr>& q, const std::vector<VectorXr>& v);
    // Assign contiguous DOF indices (3 per dynamic entity)
    void assignDofs();

    // Detect and set contacts from the current system state
    void updateContactsFromSystem();

    // Check system dirty flags and structural updates once per step and rebuild caches as needed
    void refreshState();
    void refreshCollisions();


private:
    const PhysicsSystem& m_sys;

    // Cached data
    std::vector<VectorXr> m_q;
    std::vector<VectorXr> m_v;
    std::vector<VectorXr> m_f;
    std::vector<MatrixXXr> m_Mass_blocks; // size Nb, each MatrixXXr dofVxdofV
    std::vector<MatrixXXr> m_Minv_blocks; // size Nb, each MatrixXXr dofVxdofV

    std::vector<MatrixXXr> m_W_blocks; // size Nc * 2 each MatrixXXr dofPxdofV
    std::vector<std::pair<int, int>> m_W_from_contact;
    std::vector<std::vector<int>> m_W_from_body; // size Nb, each
    std::vector<int> m_W_block_to_body; // size Nc*2, maps W-block index -> body index (-1 for static)
    std::vector<int> m_W_block_to_contact; // size Nc*2, maps W-block index -> contact id
    std::vector<collision::Contact> m_contacts;

    // Dirty flags
    bool m_contacts_dirty{true};

    // Cached sizes
    index_t m_numQ{0};
    index_t m_numV{0};

    // Rebuild helpers
    void rebuildMass_();
    void rebuildForces_();
    void rebuildW_();
    void rebuildBlocks_();
};

} // namespace cardillo::physics
