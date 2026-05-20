#include "projected_gauss_seidel.hpp"

#include <cmath>

namespace cardillo::solver {

static inline void project(VectorXr& impulse, VectorXr& dlambda, const real_t mu, const int offset, const int n) {
    const real_t impulse_normal = impulse[offset];
    impulse[offset] = std::min(impulse_normal, (real_t)0);
    dlambda[offset] += impulse[offset] - impulse_normal;

    if (n <= 1) return;

    const real_t t_norm = impulse.segment(offset + 1, n - 1).norm();

    if (t_norm <= 0) return;

    const real_t s = std::min<real_t>((real_t)1, -(mu * impulse[offset]) / t_norm);
    const Vector2r t_proj = s * impulse.segment(offset + 1, n - 1);
    dlambda.segment(offset + 1, n - 1) += t_proj - impulse.segment(offset + 1, n - 1);
    impulse.segment(offset + 1, n - 1) = t_proj;
}

VectorXr ProjectedGaussSeidelSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSetup);

    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int nSprings = m_dyn.numSprings();
    const int nDampers = m_dyn.numDampers();
    const int nContacts = m_dyn.numContactRows();
    const int Nnfc = m_dyn.numFrictionlessContacts();
    const int numLambda = nSprings + nDampers + nContacts;

    if (m_cfg.debug_pj) {
        std::cout << "PGS solve: nSprings = " << nSprings << ", nDampers = " << nDampers << ", nContacts = " << nContacts << ", numLambda = " << numLambda << std::endl;
    }

    if (numLambda == 0) return u_free;

    VectorXr mus = m_dyn.muVec();
    const auto& DinvMatrix = m_assembler.Dinv(dt, theta);
    const VectorXr rhs = m_assembler.rhs(dt, theta);
    VectorXr u_corr = VectorXr::Zero(u_free.size());

    Eigen::SparseMatrix<real_t, Eigen::RowMajor> W_row(m_assembler.W());
    Eigen::SparseMatrix<real_t, Eigen::ColMajor> W_T = W_row.transpose();  // N x M

    const VectorXr C_vec = m_assembler.C(dt, theta);

    VectorXr res = VectorXr::Zero(numLambda);
    VectorXr lambda = VectorXr::Zero(numLambda);
    VectorXr dlambda = VectorXr::Zero(numLambda);

    VectorXr projected = VectorXr::Zero(nContacts);
    VectorXr d_projected = VectorXr::Zero(nContacts);

    // Warmstart
    if (m_dyn.Lambda_g().size() != m_dyn.numSprings()) m_dyn.setLambda_g(VectorXr::Zero(m_dyn.numSprings()));
    if (m_dyn.Lambda_gamma().size() != m_dyn.numDampers()) m_dyn.setLambda_gamma(VectorXr::Zero(m_dyn.numDampers()));

    if (nSprings > 0) lambda.head(nSprings) = m_dyn.Lambda_g();
    if (nDampers > 0) lambda.tail(nDampers) = m_dyn.Lambda_gamma();

    if (m_cfg.pj_warmstart && nContacts > 0) {
        VectorXr l_contact = VectorXr::Zero(nContacts);
        WarmstartProvider::applyWarmstart(l_contact, m_dyn);
        lambda.segment(nSprings + nDampers, nContacts) = l_contact;
    }

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
            if (offset >= nSprings + nDampers) dlambda.segment(offset, blockSize) *= alpha;  // Step size scaling for contact blocks

            lambda.segment(offset, blockSize).noalias() += dlambda.segment(offset, blockSize);

            // Project lambda to satisfy constraints
            if (offset >= nSprings + nDampers) project(lambda, dlambda, mus[offset - nSprings - nDampers], offset, blockSize);

            // Update underlying u_corr via W_T
            u_corr.noalias() += m_dyn.MinvDiag().cwiseProduct(W_T.middleCols(offset, blockSize) * dlambda.segment(offset, blockSize));
            //  for (int c = 0; c < blockSize; ++c) {
            //     real_t dl = dlambda[offset + c];
            //     for (Eigen::SparseMatrix<real_t, Eigen::ColMajor>::InnerIterator it(W_T, offset + c); it; ++it) {
            //         u_corr[it.row()] += m_dyn.MinvDiag()[it.row()] * it.value() * dl;
            //     }
            // }

            // Accumulate residual norm for convergence check (use dlambda instead of res to get projected residual)
            res_norm_sq += dlambda.segment(offset, blockSize).squaredNorm() / (DinvBlocks[i].diagonal().squaredNorm());
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

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) WarmstartProvider::storeImpulse(lambda.segment(nSprings + nDampers, nContacts), m_dyn);

    return u_free - u_corr;
}

}  // namespace cardillo::solver