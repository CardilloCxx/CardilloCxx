#include "dynamics_assembler.hpp"

namespace cardillo::physics {

std::optional<index_t> DynamicsAssembler::velDofStart(entt::entity e) const {
    if (m_sys.ecs().any_of<PhysicsSystem::C_LinearVelocityIndex3>(e)) {
        const auto& vi = m_sys.ecs().get<PhysicsSystem::C_LinearVelocityIndex3>(e);
        return vi.idx.start;
    }
    return std::nullopt;
}

void DynamicsAssembler::updateContactsFromSystem() {
    m_contacts = cardillo::collision::detectAll(m_sys);
    m_contacts_dirty = true;
    m_G_dirty = true;
}

// ---------- Cached API ----------

const VectorXr& DynamicsAssembler::q() { return m_q; }

const VectorXr& DynamicsAssembler::v() { return m_v; }

const VectorXr& DynamicsAssembler::f() { return m_f; }

const Eigen::SparseMatrix<real_t>& DynamicsAssembler::M() { return m_M; }

const Eigen::SparseMatrix<real_t>& DynamicsAssembler::Minv() { return m_Minv; }

const VectorXr& DynamicsAssembler::MinvDiag() { return m_MinvDiag; }

const Eigen::SparseMatrix<real_t>& DynamicsAssembler::W() { return m_W; }

const Eigen::SparseMatrix<real_t>& DynamicsAssembler::G() { return m_G; }

void DynamicsAssembler::setContacts(std::vector<collision::Contact> contacts) {
    m_contacts = std::move(contacts);
    m_contacts_dirty = true;
    m_G_dirty = true;
}

void DynamicsAssembler::markContactsDirty() {
    m_contacts_dirty = true;
    m_G_dirty = true;
}

// ---------- Rebuild helpers ----------

void DynamicsAssembler::rebuildMass_() {
    // Build mass diagonal directly from ECS
    const index_t V = numV();
    m_Mdiag = VectorXr::Zero(V);
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_Mass, PhysicsSystem::C_LinearVelocityIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, m, vi] : view.each()) {
            (void)e;
            m_Mdiag.segment<3>(vi.idx.start).setConstant(m.m);
        }
    }
    std::vector<Eigen::Triplet<real_t>> tripsM;
    tripsM.reserve(static_cast<size_t>(V));
    for (index_t i = 0; i < V; ++i) if (m_Mdiag[i] != (real_t)0) tripsM.emplace_back(i, i, m_Mdiag[i]);
    m_M.resize(V, V);
    m_M.setFromTriplets(tripsM.begin(), tripsM.end());

    // Build inverse mass as sparse diagonal
    m_MinvDiag.resize(V);
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve(static_cast<size_t>(V));
    for (index_t i = 0; i < V; ++i) {
        real_t inv = (m_Mdiag[i] != 0) ? (1.0 / m_Mdiag[i]) : 0.0;
        if (inv != (real_t)0) trips.emplace_back(i, i, inv);
        m_MinvDiag[i] = inv;
    }
    m_Minv.resize(V, V);
    m_Minv.setFromTriplets(trips.begin(), trips.end());
    m_Minv.makeCompressed();
}

void DynamicsAssembler::rebuildForces_() {
    // Compute f from gravity and masses
    const Vector3r& g = m_sys.gravity();
    const index_t V = numV();
    m_f.resize(V);
    for (index_t i = 0; i < V; i += 3) {
        m_f[i+0] = m_Mdiag[i+0] * g[0];
        m_f[i+1] = m_Mdiag[i+1] * g[1];
        m_f[i+2] = m_Mdiag[i+2] * g[2];
    }
}

void DynamicsAssembler::loadStateFromSystem() {
    // Pack positions and velocities from ECS
    m_q = VectorXr::Zero(numQ());
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_Position3, PhysicsSystem::C_PositionIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, x, qi] : view.each()) {
            (void)e;
            m_q.segment<3>(qi.idx.start) = x.value;
        }
    }
    m_v = VectorXr::Zero(numV());
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_LinearVelocity3, PhysicsSystem::C_LinearVelocityIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, vel, vi] : view.each()) {
            (void)e;
            m_v.segment<3>(vi.idx.start) = vel.value;
        }
    }
}

void DynamicsAssembler::writeStateToSystem(const RefVectorXr& q, const RefVectorXr& v) {
    // Unpack to ECS and mark system state dirty for next rebuild
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_Position3, PhysicsSystem::C_PositionIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, x, qi] : view.each()) {
            (void)e;
            const_cast<PhysicsSystem::C_Position3&>(x).value = q.segment<3>(qi.idx.start);
        }
    }
    {
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_LinearVelocity3, PhysicsSystem::C_LinearVelocityIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, vel, vi] : view.each()) {
            (void)e;
            const_cast<PhysicsSystem::C_LinearVelocity3&>(vel).value = v.segment<3>(vi.idx.start);
        }
    }
    m_sys.markStateDirty();
}

void DynamicsAssembler::assignDofs() {
    // Append-only DOF assignment: keep existing indices stable; assign to new entities only
    const auto& reg = m_sys.ecs();
    // Determine next free index based on existing velocity indices
    index_t q_next = 0;
    {
        auto viewIdx = reg.view<PhysicsSystem::C_LinearVelocityIndex3, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, vi] : viewIdx.each()) {
            (void)e;
            q_next = std::max(q_next, vi.idx.start + (index_t)3);
        }
    }
    auto view = reg.view<PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_Position3, PhysicsSystem::C_LinearVelocity3>();
    for (auto e : view) {
        auto entity = static_cast<entt::entity>(e);
        auto& ncr = const_cast<entt::registry&>(m_sys.ecs());
        const bool hasQi = ncr.all_of<PhysicsSystem::C_PositionIndex3>(entity);
        const bool hasVi = ncr.all_of<PhysicsSystem::C_LinearVelocityIndex3>(entity);
        if (hasQi && hasVi) continue;
        auto& qi = ncr.get_or_emplace<PhysicsSystem::C_PositionIndex3>(entity).idx;
        auto& vi = ncr.get_or_emplace<PhysicsSystem::C_LinearVelocityIndex3>(entity).idx;
        qi.start = q_next;
        vi.start = q_next;
        q_next += 3;
    }
    // Update cached sizes from max assigned index
    m_numQ = q_next;
    m_numV = q_next;
}

void DynamicsAssembler::rebuildW_() {
    const int C = static_cast<int>(m_contacts.size());
    const index_t V = numV();
    std::vector<Eigen::Triplet<real_t>> wtrips;
    wtrips.reserve(static_cast<size_t>(C) * 6); // up to 3 for a and 3 for b
    for (int i = 0; i < C; ++i) {
        const auto& c = m_contacts[i];
        auto a0 = velDofStart(c.a);
        auto b0 = velDofStart(c.b);
        const Vector3r& n = c.normal; // A -> B
        if (a0) {
            wtrips.emplace_back(i, *a0 + 0, -n[0]);
            wtrips.emplace_back(i, *a0 + 1, -n[1]);
            wtrips.emplace_back(i, *a0 + 2, -n[2]);
        }
        if (b0) {
            wtrips.emplace_back(i, *b0 + 0, n[0]);
            wtrips.emplace_back(i, *b0 + 1, n[1]);
            wtrips.emplace_back(i, *b0 + 2, n[2]);
        }
    }
    m_W.resize(C, V);
    m_W.setFromTriplets(wtrips.begin(), wtrips.end());
    m_W.makeCompressed();
    m_contacts_dirty = false;
    m_G_dirty = true; // W changed, so G must be recomputed
    (void)V; // no-op
}

void DynamicsAssembler::rebuildG_() {
    // If structure changed, mass caches must be up-to-date before computing G
    if (m_contacts_dirty) rebuildW_();

    // Note: dimensions: W (C x V), Minv (V x V), W^T (V x C) -> G (C x C)
    Eigen::SparseMatrix<real_t> B = (m_W * m_Minv).eval();
    m_G = (B * m_W.transpose()).eval();
    m_G.makeCompressed();
    m_G_dirty = false;
    (void)B; // no-op
}

void DynamicsAssembler::refreshState() {

    // Detect state changes via velocity/position updates
    bool structureChanged = false;
    if (m_sys.consumeStructureDirty()) {
        structureChanged = true;
        assignDofs();
        rebuildMass_();
    }

    // Positions/velocities: if changed, reload
    if (m_sys.consumeStateDirty() || structureChanged) {
        loadStateFromSystem();
    }

    // Forces: if gravity/external changed or sizes changed, rebuild
    if (m_sys.consumeForcesDirty() || structureChanged) {
        rebuildForces_();
    }

    // Contacts and dependent matrices
    updateContactsFromSystem();
    if (m_contacts_dirty || structureChanged) rebuildW_();
    if (m_G_dirty || m_contacts_dirty || structureChanged) rebuildG_();
}

} // namespace cardillo::physics
