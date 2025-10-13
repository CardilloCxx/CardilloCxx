#include "projected_jacobi.hpp"
#include <iostream>

namespace cardillo::solver {

static inline VectorXr projected_jacobi_core(const Eigen::SparseMatrix<real_t>& G,
											 const Eigen::SparseMatrix<real_t>& W,
											 const VectorXr& b,
											 real_t alpha,
											 std::optional<RefVectorXr> p0,
											 int iters) {

    std::cout << "[PJ] W matrix: " << W << std::endl;
    std::cout << "[PJ] G matrix: " << G << std::endl;

	const int C = G.rows();
	if (C == 0) return VectorXr();
	(void)W; // silence in release

	VectorXr D(C);
	for (int i = 0; i < C; ++i) D[i] = G.coeff(i, i);
	// no-op: diag(G) stats removed for release

	VectorXr Rdiag(C);
	for (int i = 0; i < C; ++i) Rdiag[i] = (D[i] > (real_t)0) ? (alpha / D[i]) : (real_t)0;

	VectorXr p = VectorXr::Zero(C);
	if (p0.has_value() && p0->size() == C) p = p0.value();

	VectorXr y(C);
	for (int k = 0; k < iters; ++k) {
		y = (G * p);
		for (int i = 0; i < C; ++i) y[i] = Rdiag[i] * y[i] - p[i] + Rdiag[i] * b[i];
		for (int i = 0; i < C; ++i) p[i] = std::max<real_t>(0, y[i]);
	}
	// no-op: debug print removed
	return p;
}

VectorXr ProjectedJacobiSolver::iterate(std::optional<RefVectorXr> p0, int iters) {
	// Caller expected to have refreshed state
	const auto& W = m_dyn.W();
	const auto& G = m_dyn.G();
	const VectorXr& v = m_dyn.v();
	VectorXr b = W * v;
	return projected_jacobi_core(G, W, b, m_alpha, p0, iters);
}

VectorXr ProjectedJacobiSolver::iterateWithPreliminaryVelocity(const VectorXr& v_pre,
															   std::optional<RefVectorXr> p0,
															   int iters) {
	const auto& W = m_dyn.W();
	const auto& G = m_dyn.G();
	VectorXr b = W * v_pre;
    std::cout << "[PJ] b vector: " << b.transpose() << std::endl;
	return projected_jacobi_core(G, W, b, m_alpha, p0, iters);
}

} // namespace cardillo::solver