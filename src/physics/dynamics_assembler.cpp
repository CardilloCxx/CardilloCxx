
#include "dynamics_assembler.hpp"
#include "../collision/collision_manager.hpp"
#include <algorithm>

namespace cardillo::physics {

void DynamicsAssembler::updateContactsFromSystem() {
    cardillo::collision::CollisionManager mgr;
    m_contacts = mgr.detectAll(m_sys);
    m_contacts_dirty = true;
}

// ---------- Cached API (block-based) ----------

const std::vector<VectorXr>& DynamicsAssembler::q() { return m_q; }
const std::vector<VectorXr>& DynamicsAssembler::v() { return m_v; }
const std::vector<VectorXr>& DynamicsAssembler::f() { return m_f; }

const std::vector<VectorXr>& DynamicsAssembler::vBlocks() { return m_v; }
const std::vector<MatrixXXr>& DynamicsAssembler::MinvBlocks() { return m_Minv_blocks; }
const std::vector<MatrixXXr>& DynamicsAssembler::WBlocks() { return m_W_blocks; }
const std::pair<int,int>& DynamicsAssembler::WBlocksFromContact(index_t contact) { return m_W_from_contact[(size_t)contact]; }
const std::vector<int>& DynamicsAssembler::WBlocksFromBody(index_t body) { return m_W_from_body[(size_t)body]; }

void DynamicsAssembler::setContacts(std::vector<collision::Contact> contacts) {
    m_contacts = std::move(contacts);
    m_contacts_dirty = true;
}

void DynamicsAssembler::markContactsDirty() { m_contacts_dirty = true; }

// ---------- Rebuild helpers ----------

void DynamicsAssembler::rebuildMass_() {
    const int Nb = m_sys.numBodies();
    m_Mass_blocks.assign(Nb, MatrixXXr());
    m_Minv_blocks.assign(Nb, MatrixXXr());
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        (void)e;
        const int b = bi.b;
        if (b < 0 || b >= Nb) continue;
        m_Mass_blocks[(size_t)b] = m_sys.getMass(e);
        m_Minv_blocks[(size_t)b] = m_sys.getMassInverse(e);
    }
}

void DynamicsAssembler::rebuildForces_() {
    const int Nb = m_sys.numBodies();
    m_f.assign(Nb, VectorXr());
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < Nb) m_f[(size_t)b] = m_sys.getForce(e);
    }
}

void DynamicsAssembler::loadStateFromSystem() {
    const int Nb = m_sys.numBodies();
    m_q.assign(Nb, VectorXr());
    m_v.assign(Nb, VectorXr());
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b >= 0 && b < Nb) m_q[(size_t)b] = m_sys.getPosition(e);
        }
    }
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b >= 0 && b < Nb) m_v[(size_t)b] = m_sys.getVelocity(e);
        }
    }
}

void DynamicsAssembler::writeStateToSystem(const std::vector<VectorXr>& q, const std::vector<VectorXr>& v) {
    
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_LinearVelocity3, PhysicsSystem::C_Position3, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi, vel, x] : view.each()) {
        (void)e;
        const int b = bi.b;
        if (b >= 0 && b < (int)q.size())
            const_cast<PhysicsSystem::C_Position3&>(x).value = q[(size_t)b];
        if (b >= 0 && b < (int)v.size())
            const_cast<PhysicsSystem::C_LinearVelocity3&>(vel).value = v[(size_t)b];
    }
    m_sys.markStateDirty();
}

void DynamicsAssembler::assignDofs() {
    auto& reg = const_cast<entt::registry&>(m_sys.ecs());
    // Assign consecutive body indices to dynamic entities
    int nextBody = 0;
    auto view = reg.view<PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_Position3, PhysicsSystem::C_LinearVelocity3>();
    for (auto e : view) {
        entt::entity ent = static_cast<entt::entity>(e);
        if (!reg.any_of<PhysicsSystem::C_BodyIndex>(ent)) {
            reg.emplace<PhysicsSystem::C_BodyIndex>(ent, PhysicsSystem::C_BodyIndex{nextBody});
            ++nextBody;
        }
    }
    // Update cached sizes; numQ/numV are now sums over per-entity DOFs (queried as needed)
    m_numQ = 0; m_numV = 0;
    for (auto e : view) {
        m_numQ += (index_t)m_sys.getPosition(static_cast<entt::entity>(e)).size();
        m_numV += (index_t)m_sys.getVelocity(static_cast<entt::entity>(e)).size();
    }
}

void DynamicsAssembler::rebuildW_() {
    const int C = (int)m_contacts.size();
    const int Nb = m_sys.numBodies();
    m_W_blocks.clear(); m_W_blocks.reserve((size_t)C * 2);
    m_W_from_contact.clear(); m_W_from_contact.reserve((size_t)C);
    m_W_from_body.assign((size_t)Nb, {});
    m_W_block_to_body.clear(); m_W_block_to_body.reserve((size_t)C * 2);
    m_W_block_to_contact.clear(); m_W_block_to_contact.reserve((size_t)C * 2);

    const auto& reg = m_sys.ecs();
    for (int i = 0; i < C; ++i) {
        const auto& c = m_contacts[i];
        bool aDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.a);
        bool bDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.b);

        int idxA = -1;
        if (aDyn) {
            idxA = (int)m_W_blocks.size();
            m_W_blocks.push_back(c.wA);
            m_W_block_to_body.push_back(reg.get<PhysicsSystem::C_BodyIndex>(c.a).b);
            m_W_block_to_contact.push_back(i);
            int b = reg.get<PhysicsSystem::C_BodyIndex>(c.a).b;
            if (b >= 0 && b < Nb) m_W_from_body[(size_t)b].push_back(idxA);
        }

        int idxB = -1;
        if (bDyn) {
            idxB = (int)m_W_blocks.size();
            m_W_blocks.push_back(c.wB);
            m_W_block_to_body.push_back(reg.get<PhysicsSystem::C_BodyIndex>(c.b).b);
            m_W_block_to_contact.push_back(i);
            int b = reg.get<PhysicsSystem::C_BodyIndex>(c.b).b;
            if (b >= 0 && b < Nb) m_W_from_body[(size_t)b].push_back(idxB);
        }

        m_W_from_contact.emplace_back(idxA, idxB);
    }
    m_contacts_dirty = false;
}

void DynamicsAssembler::rebuildBlocks_() {
    const int Nb = m_sys.numBodies();
    if ((int)m_q.size() != Nb) m_q.assign(Nb, VectorXr());
    if ((int)m_v.size() != Nb) m_v.assign(Nb, VectorXr());
    if ((int)m_f.size() != Nb) m_f.assign(Nb, VectorXr());
}

void DynamicsAssembler::refreshState() {
    bool structureChanged = false;
    if (m_sys.consumeStructureDirty()) {
        structureChanged = true;
        assignDofs();
        rebuildBlocks_();
        rebuildMass_();
    }

    if (m_sys.consumeStateDirty() || structureChanged) {
        loadStateFromSystem();
    }

    if (m_sys.consumeForcesDirty() || structureChanged) {
        rebuildForces_();
    }
}

void DynamicsAssembler::refreshCollisions() {
    updateContactsFromSystem();
    rebuildW_();
}


} // namespace cardillo::physics
