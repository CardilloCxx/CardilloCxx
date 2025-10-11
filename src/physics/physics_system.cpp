#include "physics_system.hpp"

namespace cardillo {

// Construction and configuration
PhysicsSystem::PhysicsSystem() {
    PetscPrintf(PETSC_COMM_WORLD, "Hello from PETSc + Eigen!\n");
    m_gravity = Vector3r(0, 0, -9.81);
    m_q_dofs = m_v_dofs = 0;
}

void PhysicsSystem::setGravity(const Vector3r& g) { m_gravity = g; }

// Entity creation
index_t PhysicsSystem::addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0) {
    auto e = m_reg.create();
    m_reg.emplace<C_PhysicsObject>(e);
    m_reg.emplace<C_PointMassTag>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_PointVisualTag>(e);
    m_reg.emplace<C_Mass>(e, C_Mass{mass});
    m_reg.emplace<C_Position3>(e, C_Position3{x0});
    m_reg.emplace<C_LinearVelocity3>(e, C_LinearVelocity3{v0});
    assignDofs_();
    m_mass_dirty = true;
    return static_cast<index_t>(entt::to_integral(e));
}

entt::entity PhysicsSystem::createRigidVisualEntity_(const Vector3r& center) {
    auto e = m_reg.create();
    m_reg.emplace<C_RigidBodyTag>(e);
    m_reg.emplace<C_VisualObject>(e);
    m_reg.emplace<C_Position3>(e, C_Position3{center});
    return e;
}

index_t PhysicsSystem::addRigidBody(const Plane& p) {
    auto e = createRigidVisualEntity_(p.center);
    m_reg.emplace<C_PlaneVisualTag>(e);
    m_reg.emplace<C_Plane>(e, C_Plane{p.normal, p.up, p.sizeX, p.sizeY});
    return static_cast<index_t>(entt::to_integral(e));
}

index_t PhysicsSystem::addRigidBody(const Cube& c) {
    auto e = createRigidVisualEntity_(c.center);
    m_reg.emplace<C_CubeVisualTag>(e);
    m_reg.emplace<C_Cube>(e, C_Cube{c.halfExtents});
    return static_cast<index_t>(entt::to_integral(e));
}

// DOF management
void PhysicsSystem::assignDofs_() {
    index_t q_next = 0;
    // assign or update indices on all dynamic bodies
    auto view = m_reg.view<C_PhysicsObject, C_Position3, C_LinearVelocity3>();
    for (auto e : view) {
        auto& qi = m_reg.get_or_emplace<C_PositionIndex3>(e).idx;
        auto& vi = m_reg.get_or_emplace<C_LinearVelocityIndex3>(e).idx;
        qi.start = q_next;
        vi.start = q_next;
        q_next += 3;
    }
    m_q_dofs = m_v_dofs = q_next;
    m_mass_dirty = true;
}

// Assembly
Eigen::SparseMatrix<real_t> PhysicsSystem::assembleMassMatrix() const {
    std::vector<Eigen::Triplet<real_t>> trips;
    size_t count = 0;
    auto view = m_reg.view<C_Mass, C_LinearVelocityIndex3, C_PhysicsObject>();
    for (auto [e, m, vi] : view.each()) {
        (void)e; (void)m; (void)vi;
        count += 3;
    }
    trips.reserve(count);
    for (auto [e, m, vi] : view.each()) {
        (void)e;
        for (int i = 0; i < 3; ++i) trips.emplace_back(vi.idx.start + i, vi.idx.start + i, m.m);
    }
    Eigen::SparseMatrix<real_t> M(numV(), numV());
    M.setFromTriplets(trips.begin(), trips.end());
    return M;
}

VectorXr PhysicsSystem::assembleForceVector() const {
    VectorXr f = VectorXr::Zero(numV());
    auto view = m_reg.view<C_Mass, C_LinearVelocityIndex3, C_PhysicsObject>();
    for (auto [e, m, vi] : view.each()) {
        (void)e;
        f.segment<3>(vi.idx.start) = m.m * m_gravity;
    }
    return f;
}

// State pack/unpack
VectorXr PhysicsSystem::packQ() const {
    VectorXr q = VectorXr::Zero(numQ());
    auto view = m_reg.view<C_Position3, C_PositionIndex3, C_PhysicsObject>();
    for (auto [e, x, qi] : view.each()) {
        (void)e;
        q.segment<3>(qi.idx.start) = x.value;
    }
    return q;
}

VectorXr PhysicsSystem::packV() const {
    VectorXr v = VectorXr::Zero(numV());
    auto view = m_reg.view<C_LinearVelocity3, C_LinearVelocityIndex3, C_PhysicsObject>();
    for (auto [e, vel, vi] : view.each()) {
        (void)e;
        v.segment<3>(vi.idx.start) = vel.value;
    }
    return v;
}

void PhysicsSystem::unpackQ(const RefVectorXr& q) {
    auto view = m_reg.view<C_Position3, C_PositionIndex3, C_PhysicsObject>();
    for (auto [e, x, qi] : view.each()) {
        (void)e;
        x.value = q.segment<3>(qi.idx.start);
    }
}

void PhysicsSystem::unpackV(const RefVectorXr& v) {
    auto view = m_reg.view<C_LinearVelocity3, C_LinearVelocityIndex3, C_PhysicsObject>();
    for (auto [e, vel, vi] : view.each()) {
        (void)e;
        vel.value = v.segment<3>(vi.idx.start);
    }
}

// Cached mass diagonal
const VectorXr& PhysicsSystem::massDiagonal() const {
    if (!m_mass_dirty && m_Mdiag.size() == numV()) return m_Mdiag;
    m_Mdiag = VectorXr::Zero(numV());
    auto view = m_reg.view<C_Mass, C_LinearVelocityIndex3, C_PhysicsObject>();
    for (auto [e, m, vi] : view.each()) {
        (void)e;
        m_Mdiag.segment<3>(vi.idx.start).setConstant(m.m);
    }
    m_mass_dirty = false;
    return m_Mdiag;
}

} // namespace cardillo
