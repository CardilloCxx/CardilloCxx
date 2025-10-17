
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
    // Point masses
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_PointMassTag>();
        for (auto [e, bi] : view.each()) {
            (void)e;
            const int b = bi.b; if (b < 0 || b >= Nb) continue;
            m_Mass_blocks[(size_t)b] = m_sys.getMass(e);
            m_Minv_blocks[(size_t)b] = m_sys.getMassInverse(e);
        }
    }
    // Rigid bodies
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_RigidBodyTag>();
        for (auto [e, bi] : view.each()) {
            (void)e;
            const int b = bi.b; if (b < 0 || b >= Nb) continue;
            m_Mass_blocks[(size_t)b] = m_sys.getMass(e);
            m_Minv_blocks[(size_t)b] = m_sys.getMassInverse(e);
        }
    }
}

void DynamicsAssembler::rebuildForces_() {
    const int Nb = m_sys.numBodies();
    m_f.assign(Nb, VectorXr());
    const auto& reg = m_sys.ecs();
    // Point masses
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_PointMassTag>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b; if (b >= 0 && b < Nb) m_f[(size_t)b] = m_sys.getForce(e);
        }
    }
    // Rigid bodies
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_RigidBodyTag>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b; if (b >= 0 && b < Nb) m_f[(size_t)b] = m_sys.getForce(e);
        }
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
    // Point masses: q size 3, v size 3
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_PointMassTag, PhysicsSystem::C_LinearVelocity3, PhysicsSystem::C_Position3>();
        for (auto [e, bi, vel, x] : view.each()) {
            (void)e;
            const int b = bi.b;
            if (b >= 0 && b < (int)q.size() && q[(size_t)b].size() >= 3) {
                const VectorXr& qb = q[(size_t)b];
                const_cast<PhysicsSystem::C_Position3&>(x).value = Vector3r(qb[0], qb[1], qb[2]);
            }
            if (b >= 0 && b < (int)v.size() && v[(size_t)b].size() >= 3) {
                const VectorXr& vb = v[(size_t)b];
                const_cast<PhysicsSystem::C_LinearVelocity3&>(vel).value = Vector3r(vb[0], vb[1], vb[2]);
            }
        }
    }
    // Rigid bodies: q size 7 ([x(3), w, px, py, pz]), v size 6 ([v(3), omega_body(3)])
    {
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_RigidBodyTag,
                              PhysicsSystem::C_LinearVelocity3, PhysicsSystem::C_AngularVelocity3,
                              PhysicsSystem::C_Position3, PhysicsSystem::C_Orientation>();
        for (auto [e, bi, vlin, vang, x, ori] : view.each()) {
            (void)e;
            const int b = bi.b;
            if (b >= 0 && b < (int)q.size() && q[(size_t)b].size() >= 7) {
                const VectorXr& qb = q[(size_t)b];
                Vector3r xn(qb[0], qb[1], qb[2]);
                Quaternion4r qn(qb[3], qb[4], qb[5], qb[6]);
                qn.normalize();
                const_cast<PhysicsSystem::C_Position3&>(x).value = xn;
                const_cast<PhysicsSystem::C_Orientation&>(ori).q = qn;
            }
            if (b >= 0 && b < (int)v.size() && v[(size_t)b].size() >= 6) {
                const VectorXr& vb = v[(size_t)b];
                const_cast<PhysicsSystem::C_LinearVelocity3&>(vlin).value = Vector3r(vb[0], vb[1], vb[2]);
                const_cast<PhysicsSystem::C_AngularVelocity3&>(vang).value = Vector3r(vb[3], vb[4], vb[5]);
            }
        }
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
    m_contact_index_orig.clear(); m_contact_index_orig.reserve((size_t)C);

    const auto& reg = m_sys.ecs();
    int dynContactId = 0; // index into dynamic contacts (those with at least one dynamic body)
    for (int i = 0; i < C; ++i) {
        const auto& c = m_contacts[i];
        bool aDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.a);
        bool bDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.b);

        // Skip pure static-static contacts for W assembly; they don't contribute to impulses
        if (!aDyn && !bDyn) {
            continue;
        }

        int idxA = -1;
        if (aDyn) {
            idxA = (int)m_W_blocks.size();
            // Extend W row for rigid body A: [ -n, r_A x (-n) ] where r_A = p - x_A (world), torque term expressed in body frame.
            // Justification: y = n·(v_B + ω_B×r_B − v_A − ω_A×r_A) =
            //   (-n)·v_A  + (−(r_A×n))·ω_A   +   (+n)·v_B  + ((r_B×n))·ω_B
            MatrixXXr row = c.wA;
            if (reg.any_of<PhysicsSystem::C_RigidBodyTag>(c.a)) {
                // Build 1x6 row
                MatrixXXr w(1,6); w.setZero();
                // translational part
                w(0,0) = row(0,0); if (row.cols() > 1) w(0,1) = row(0,1); if (row.cols() > 2) w(0,2) = row(0,2);
                // torque part for A using body-space values from contact: r_body × (−n_body)
                Vector3r t_body = c.pointA_body.cross(-c.normalA_body);
                w(0,3) = t_body.x(); w(0,4) = t_body.y(); w(0,5) = t_body.z();
                row = w;
            }
            m_W_blocks.push_back(row);
            m_W_block_to_body.push_back(reg.get<PhysicsSystem::C_BodyIndex>(c.a).b);
            m_W_block_to_contact.push_back(dynContactId);
            int b = reg.get<PhysicsSystem::C_BodyIndex>(c.a).b;
            if (b >= 0 && b < Nb) m_W_from_body[(size_t)b].push_back(idxA);
        }

        int idxB = -1;
        if (bDyn) {
            idxB = (int)m_W_blocks.size();
            MatrixXXr row = c.wB;
            if (reg.any_of<PhysicsSystem::C_RigidBodyTag>(c.b)) {
                MatrixXXr w(1,6); w.setZero();
                w(0,0) = row(0,0); if (row.cols() > 1) w(0,1) = row(0,1); if (row.cols() > 2) w(0,2) = row(0,2);
                // torque part for B using body-space values from contact: r_body × (n_body)
                Vector3r t_body = c.pointB_body.cross(c.normalB_body);
                w(0,3) = t_body.x(); w(0,4) = t_body.y(); w(0,5) = t_body.z();
                row = w;
            }
            m_W_blocks.push_back(row);
            m_W_block_to_body.push_back(reg.get<PhysicsSystem::C_BodyIndex>(c.b).b);
            m_W_block_to_contact.push_back(dynContactId);
            int b = reg.get<PhysicsSystem::C_BodyIndex>(c.b).b;
            if (b >= 0 && b < Nb) m_W_from_body[(size_t)b].push_back(idxB);
        }

        m_W_from_contact.emplace_back(idxA, idxB);
        m_contact_index_orig.push_back(i);
        ++dynContactId;
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
