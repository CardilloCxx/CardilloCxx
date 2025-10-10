#include "world.hpp"

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
	// Simple contiguous layout: each mass contributes 3 DoFs for q and v
	index_t q_next = 0;
	for (auto& pm : m_masses) {
		pm.q_start = q_next;
		pm.v_start = q_next; // identical for translational-only
		q_next += 3;
	}
	m_q_dofs = m_v_dofs = static_cast<size_t>(q_next);
}

Eigen::SparseMatrix<real_t> PhysicsSystem::assembleMassMatrix() const
{
	std::vector<Eigen::Triplet<real_t>> trips;
	trips.reserve(3 * m_masses.size());
	for (const auto& pm : m_masses) {
		for (int i = 0; i < 3; ++i) {
			trips.emplace_back(pm.v_start + i, pm.v_start + i, pm.m);
		}
	}
	Eigen::SparseMatrix<real_t> M(numV(), numV());
	M.setFromTriplets(trips.begin(), trips.end());
	return M;
}

VectorXr PhysicsSystem::assembleForceVector() const
{
	VectorXr f = VectorXr::Zero(numV());
	for (const auto& pm : m_masses) {
		f.segment<3>(pm.v_start) = pm.m * m_gravity;
	}
	return f;
}

VectorXr PhysicsSystem::packQ() const
{
	VectorXr q(numQ());
	for (const auto& pm : m_masses) {
		q.segment<3>(pm.q_start) = pm.x;
	}
	return q;
}

VectorXr PhysicsSystem::packV() const
{
	VectorXr v(numV());
	for (const auto& pm : m_masses) {
		v.segment<3>(pm.v_start) = pm.v;
	}
	return v;
}

void PhysicsSystem::unpackQ(const RefVectorXr& q)
{
	for (auto& pm : m_masses) {
		pm.x = q.segment<3>(pm.q_start);
	}
}

void PhysicsSystem::unpackV(const RefVectorXr& v)
{
	for (auto& pm : m_masses) {
		pm.v = v.segment<3>(pm.v_start);
	}
}

void PhysicsSystem::stepMidpoint(real_t dt)
{
	// Pack state
	VectorXr qn = packQ();
	VectorXr vn = packV();

	// Mass matrix and its inverse solve prepared once (diagonal so easy)
	// Compute acceleration at midpoint using a predictor for midpoint state:
	// predictor: v_mid ~ v_n + 0.5*dt * a_n; q_mid ~ q_n + 0.5*dt * v_n
	Eigen::SparseMatrix<real_t> M = assembleMassMatrix();

	// a_n = M^{-1} f_n
	VectorXr fn = assembleForceVector();
	// Solve M * a = f (M diagonal -> direct division)
	VectorXr an(numV());
	an.setZero();
	// Extract diagonal once
	VectorXr Mdiag = VectorXr::Zero(numV());
	for (int k = 0; k < M.outerSize(); ++k) {
		for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(M, k); it; ++it) {
			if (it.row() == it.col()) Mdiag[it.row()] = it.value();
		}
	}
	for (index_t i = 0; i < numV(); ++i) an[i] = fn[i] / Mdiag[i];

	VectorXr v_mid = vn + (real_t)0.5 * dt * an;
	VectorXr q_mid = qn + (real_t)0.5 * dt * vn;

	// Evaluate forces at midpoint (gravity independent of state here, but keep API)
	VectorXr f_mid = assembleForceVector();

	// a_mid = M^{-1} f_mid
	VectorXr a_mid(numV());
	for (index_t i = 0; i < numV(); ++i) a_mid[i] = f_mid[i] / Mdiag[i];

	// Update to n+1
	VectorXr vnp1 = vn + dt * a_mid;
	VectorXr qnp1 = qn + dt * v_mid;

	// Scatter back to entities
	unpackQ(qnp1);
	unpackV(vnp1);
}

} // namespace cardillo