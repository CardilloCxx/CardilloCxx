#include "projected_gauss_seidel.hpp"

#include <cmath>

namespace cardillo::solver {

VectorXr ProjectedGaussSeidelSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSetup);

    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int numSprings = m_dyn.numSprings();
    const int numDampers = m_dyn.numDampers();
    if (numSprings + numDampers == 0) return u_free;

    const auto& DinvMatrix = m_assembler.DinvDiag(dt, theta);
    const VectorXr rhs = m_assembler.rhs(dt, theta);
    VectorXr u_corr = VectorXr::Zero(u_free.size());

    Eigen::SparseMatrix<real_t, Eigen::RowMajor> W_row(m_assembler.W());
    Eigen::SparseMatrix<real_t, Eigen::ColMajor> W_T = W_row.transpose();  // N x M

    const VectorXr C_vec = m_assembler.C(dt, theta);

    VectorXr res = VectorXr::Zero(numSprings + numDampers);
    VectorXr lambda = VectorXr::Zero(numSprings + numDampers);
    VectorXr dlambda = VectorXr::Zero(numSprings + numDampers);

    // Warmstart
    if (m_dyn.Lambda_g().size() != numSprings) m_dyn.setLambda_g(VectorXr::Zero(numSprings));
    if (m_dyn.Lambda_gamma().size() != numDampers) m_dyn.setLambda_gamma(VectorXr::Zero(numDampers));

    if (numSprings > 0) lambda.head(numSprings) = m_dyn.Lambda_g();
    if (numDampers > 0) lambda.tail(numDampers) = m_dyn.Lambda_gamma();

    auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidel);

    const real_t alpha = m_cfg.pj_alpha;

    // Initial u_corr from warmstarted lambda
    u_corr = m_dyn.MinvDiag().cwiseProduct(W_T * lambda);

    const auto& DinvBlocks = DinvMatrix.blocks();
    const int nBlocks = (int)DinvBlocks.size();

    std::vector<int> offsets(nBlocks);
    for (int i = 0, o = 0; i < nBlocks; ++i) {
        offsets[i] = o;
        o += DinvBlocks[i].rows();
    }

    for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSweep);

        real_t res_norm_sq = 0;
        const bool backward = iter & 1;

        for (int bi = 0; bi < nBlocks; ++bi) {
            int i = backward ? (nBlocks - 1 - bi) : bi;
            int blockSize = DinvBlocks[i].rows();
            int offset = offsets[i];

            // Compute residual for this block
            res.segment(offset, blockSize).noalias() = rhs.segment(offset, blockSize);
            res.segment(offset, blockSize).noalias() -= W_row.middleRows(offset, blockSize) * u_corr;
            res.segment(offset, blockSize).noalias() -= C_vec.segment(offset, blockSize).cwiseProduct(lambda.segment(offset, blockSize));

            // Compute lambda update for this block
            dlambda.segment(offset, blockSize).noalias() = (DinvBlocks[i] * res.segment(offset, blockSize)) * m_cfg.pj_relaxation;
            lambda.segment(offset, blockSize).noalias() += dlambda.segment(offset, blockSize);

            // Update underlying u_corr via W_T
            u_corr.noalias() += m_dyn.MinvDiag().cwiseProduct(W_T.middleCols(offset, blockSize) * dlambda.segment(offset, blockSize));

            // Accumulate residual norm for convergence check
            res_norm_sq += res.segment(offset, blockSize).squaredNorm();
        }

        real_t res_norm = std::sqrt(res_norm_sq);

        if (std::isnan(res_norm) || std::isinf(res_norm)) {
            std::cerr << "[ProjectedGaussSeidel] Divergence detected! res_norm = " << res_norm << "\n";
            throw std::runtime_error("Projected Gauss-Seidel diverged");
        }

        if (m_cfg.debug_pj && iter % 1000 == 0) std::cout << "PGS iter " << iter << ", residual norm: " << res_norm << std::endl;

        if (res_norm < m_cfg.pj_tol_abs) break;
        m_last_iters = iter + 1;
    }

    if (numSprings > 0) m_dyn.setLambda_g(lambda.head(numSprings));
    if (numDampers > 0) m_dyn.setLambda_gamma(lambda.tail(numDampers));

    return u_free - u_corr;
}

}  // namespace cardillo::solver