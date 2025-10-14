#include "projected_jacobi.hpp"
#include <iostream>

namespace cardillo::solver {

static inline VectorXr projected_jacobi_core(const Eigen::SparseMatrix<real_t>& G,
											 const Eigen::SparseMatrix<real_t>& W,
											 const VectorXr& b,
											 real_t alpha,
											 std::optional<RefVectorXr> p0,
											 real_t tol) {

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
    VectorXr p_prev(C); 

    real_t err = std::numeric_limits<real_t>::max();
    index_t iteration = 0;
	while (err > tol && iteration < 100000) {
		y = (G * p);
        p_prev = p;
        for (int i = 0; i < C; ++i) p[i] = -std::min<real_t>(0, Rdiag[i] * y[i] - p[i] + Rdiag[i] * b[i]);

        err = (p_prev - p).norm();
        ++iteration;
        if (iteration % 1000 == 0) {
            std::cout << "[PJ] iteration " << iteration << ", err = " << err << "\n";
        }
	}

    // std::cout << "[PJ] final iteration " << iteration << ", err = " << err << "\n";
	return p;
}

VectorXr ProjectedJacobiSolver::iterate(std::optional<RefVectorXr> p0, real_t tol) {
	// Caller expected to have refreshed state
	const auto& W = m_dyn.W();
	const auto& G = m_dyn.G();
	const VectorXr& v = m_dyn.v();
	VectorXr b = W * v;
	return projected_jacobi_core(G, W, b, m_alpha, p0, tol);
}

VectorXr ProjectedJacobiSolver::iterateWithPreliminaryVelocity(const VectorXr& v_pre,
															   std::optional<RefVectorXr> p0,
															   real_t tol) {
	const auto& W = m_dyn.W();
	const auto& G = m_dyn.G();
	VectorXr b = W * v_pre;
	return projected_jacobi_core(G, W, b, m_alpha, p0, tol);
}

} // namespace cardillo::solver