#include "physics_system.hpp"

namespace cardillo {

PhysicsSystem::PhysicsSystem()
{
    m_gravity = Vector3r(0, 0, -9.81);
    m_q_dofs = m_v_dofs = 0;
}

void PhysicsSystem::setGravity(const Vector3r& g)
{
    m_gravity = g;
}

index_t PhysicsSystem::addPointMass(real_t mass, const Vector3r& x0, const Vector3r& v0)
{
    PointMass pm;
    pm.m = mass;
    pm.x = x0;
    pm.v = v0;
    m_masses.push_back(pm);
    assignDofs_();
    return static_cast<index_t>(m_masses.size() - 1);
}

void PhysicsSystem::assignDofs_()
{
    index_t q_next = 0;
    for (auto& pm : m_masses) {
        pm.q_start = q_next;
        pm.v_start = q_next;
        q_next += 3;
    }
    m_q_dofs = m_v_dofs = static_cast<size_t>(q_next);
}

Eigen::SparseMatrix<real_t> PhysicsSystem::assembleMassMatrix() const
{
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve(3 * m_masses.size());
    for (const auto& pm : m_masses) {
        for (int i = 0; i < 3; ++i) trips.emplace_back(pm.v_start + i, pm.v_start + i, pm.m);
    }
    Eigen::SparseMatrix<real_t> M(numV(), numV());
    M.setFromTriplets(trips.begin(), trips.end());
    return M;
}

VectorXr PhysicsSystem::assembleForceVector() const
{
    VectorXr f = VectorXr::Zero(numV());
    for (const auto& pm : m_masses) f.segment<3>(pm.v_start) = pm.m * m_gravity;
    return f;
}

VectorXr PhysicsSystem::packQ() const
{
    VectorXr q(numQ());
    for (const auto& pm : m_masses) q.segment<3>(pm.q_start) = pm.x;
    return q;
}

VectorXr PhysicsSystem::packV() const
{
    VectorXr v(numV());
    for (const auto& pm : m_masses) v.segment<3>(pm.v_start) = pm.v;
    return v;
}

void PhysicsSystem::unpackQ(const RefVectorXr& q)
{
    for (auto& pm : m_masses) pm.x = q.segment<3>(pm.q_start);
}

void PhysicsSystem::unpackV(const RefVectorXr& v)
{
    for (auto& pm : m_masses) pm.v = v.segment<3>(pm.v_start);
}

} // namespace cardillo
