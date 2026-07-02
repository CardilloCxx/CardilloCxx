#include "projected_gauss_seidel.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace cardillo::solver {

static inline void project(VectorXr& impulse, VectorXr& dlambda, const real_t mu, const int offset, const int n) {
    const real_t impulse_normal = impulse[offset];
    impulse[offset] = std::min(impulse_normal, (real_t)0);
    dlambda[offset] += impulse[offset] - impulse_normal;

    if (n <= 1) return;

    const VectorXr tangential = impulse.segment(offset + 1, n - 1);
    const real_t t_norm = tangential.norm();
    if (t_norm <= (real_t)0) return;

    const real_t s = std::min<real_t>((real_t)1, -(mu * impulse[offset]) / t_norm);
    const VectorXr t_proj = s * tangential;
    dlambda.segment(offset + 1, n - 1) += t_proj - tangential;
    impulse.segment(offset + 1, n - 1) = t_proj;
}

static inline void project_all(VectorXr& impulse, VectorXr& dlambda, const VectorXr& mu, const int nSprings, const int nDampers, int Nnfc, int Nfc) {
    const int nContacts = Nnfc + Nfc;
    const int nDyn = nSprings + nDampers;

    if (nContacts == 0) return;

    VectorXr impulse_normal = impulse.segment(nDyn, nContacts);

    impulse.segment(nDyn, Nnfc) = impulse.segment(nDyn, Nnfc).cwiseMin((real_t)0);
    dlambda.segment(nDyn, Nnfc) += impulse.segment(nDyn, Nnfc) - impulse_normal.head(Nnfc);

    for (int i = 0; i < Nfc; i += 3) {
        int idx = nDyn + Nnfc + i;
        const real_t impulse_n = impulse[idx];
        impulse[idx] = std::min(impulse_n, (real_t)0);
        dlambda[idx] += impulse[idx] - impulse_n;

        const VectorXr tangential = impulse.segment(idx + 1, 2);
        const real_t t_norm = tangential.norm();
        if (t_norm <= (real_t)0) continue;

        const real_t s = std::min<real_t>((real_t)1, -(mu[nDyn + Nnfc + i] * impulse[idx]) / t_norm);
        const VectorXr t_proj = s * tangential;
        dlambda.segment(idx + 1, 2) += t_proj - tangential;
        impulse.segment(idx + 1, 2) = t_proj;
    }
}

VectorXr ProjectedGaussSeidelSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSetup);

    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int nSprings = m_dyn.numSprings();
    const int nDampers = m_dyn.numDampers();
    const int nContacts = m_dyn.numContactRows();
    const int Nnfc = m_dyn.numFrictionlessContacts();
    const int Nfc = nContacts - Nnfc;
    const int nV = u_free.size();
    const int numLambda = nSprings + nDampers + nContacts;

    if (m_cfg.debug_pj) {
        std::cout << "PGS solve: nSprings = " << nSprings << ", nDampers = " << nDampers << ", nContacts = " << nContacts << ", numLambda = " << numLambda << std::endl;
    }

    if (numLambda == 0) return u_free;

    const VectorXr mus = m_dyn.muVec();
    const auto& DinvMatrix = m_assembler.DinvDiag(dt, theta);
    const VectorXr rhs = m_assembler.rhs(dt, theta);
    VectorXr u_corr = VectorXr::Zero(nV);

    Eigen::SparseMatrix<real_t, Eigen::RowMajor> W_row(m_assembler.W());
    Eigen::SparseMatrix<real_t, Eigen::ColMajor> MiW_T = m_dyn.MinvDiag().asDiagonal() * W_row.transpose();

    const VectorXr C_vec = m_assembler.C(dt, theta);

    VectorXr res = VectorXr::Zero(numLambda);
    VectorXr lambda = VectorXr::Zero(numLambda);
    VectorXr dlambda = VectorXr::Zero(numLambda);

    // Warmstart
    if (m_dyn.Lambda_g().size() != m_dyn.numSprings()) m_dyn.setLambda_g(VectorXr::Zero(m_dyn.numSprings()));
    if (m_dyn.Lambda_gamma().size() != m_dyn.numDampers()) m_dyn.setLambda_gamma(VectorXr::Zero(m_dyn.numDampers()));
    if (nSprings > 0) lambda.head(nSprings) = m_dyn.Lambda_g();
    if (nDampers > 0) lambda.segment(nSprings, nDampers) = m_dyn.Lambda_gamma();

    if (m_cfg.pj_warmstart && nContacts > 0) {
        VectorXr l_contact = VectorXr::Zero(nContacts);
        // PGS's internal convention is normal impulse <= 0 (see project_all() above); invert at
        // this boundary to/from the canonical non-negative storage convention (see warmstart.hpp).
        WarmstartProvider::applyWarmstart(l_contact, m_dyn, /*invertNormalSign=*/true);
        lambda.segment(nSprings + nDampers, nContacts) = l_contact;
    }

    auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidel);

    const real_t alpha = m_cfg.pj_alpha;

    // Initial u_corr from warmstarted lambda
    u_corr = MiW_T * lambda;

    const auto& DinvBlocks = DinvMatrix.blocks();
    const int nBlocks = (int)DinvBlocks.size();

    std::vector<int> offsets(nBlocks);
    for (int i = 0, o = 0; i < nBlocks; ++i) {
        offsets[i] = o;
        o += DinvBlocks[i].rows();
    }

    VectorXr res_diff(nV);
    VectorXr res_q(nV);
    VectorXr res_v0(nV);
    VectorXr res_v1(nV);

    auto residual = [&](const VectorXr& u_corr_k, const VectorXr& u_corr_k1) -> real_t {
        const int Nv = (int)u_corr_k.size();
        res_diff.noalias() = u_corr_k1 - u_corr_k;
        res_v0.noalias() = u_free - u_corr_k;
        res_v1.noalias() = u_free - u_corr_k1;
        res_q.array() = res_v0.cwiseAbs().cwiseMax(res_v1.cwiseAbs()).array() * m_cfg.pj_tol_rel + m_cfg.pj_tol_abs;
        return (real_t)(1.0 / std::sqrt((double)Nv) * res_diff.cwiseQuotient(res_q).norm());
    };

    auto check_residual = [&](const VectorXr& u_k, const VectorXr& u_k1, int iter) -> real_t {
        const real_t res_norm = residual(u_k1, u_k);
        if (std::isnan(res_norm) || std::isinf(res_norm)) {
            std::cerr << "[ProjectedGaussSeidel] Divergence detected after " << iter << " iterations. Residual norm: " << res_norm << std::endl;
            throw std::runtime_error("Projected Gauss-Seidel diverged");
        }
        if (m_cfg.debug_pj && iter % 1000 == 0) std::cout << "PGS iter " << iter << ", residual norm: " << res_norm << std::endl;
        return res_norm;
    };

    // Single sweep callable used by both standard and Nesterov loops
    auto pgs_sweep = [&](int iter) -> void {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSweep);
        // real_t res_norm_sq = 0;
        const bool backward = false;

        for (int bi = 0; bi < nBlocks; ++bi) {
            int i = backward ? (nBlocks - 1 - bi) : bi;
            int blockSize = DinvBlocks[i].rows();
            int offset = offsets[i];

            res.segment(offset, blockSize).noalias() = rhs.segment(offset, blockSize);
            res.segment(offset, blockSize).noalias() -= W_row.middleRows(offset, blockSize) * u_corr;
            res.segment(offset, blockSize).noalias() -= C_vec.segment(offset, blockSize).cwiseProduct(lambda.segment(offset, blockSize));

            dlambda.segment(offset, blockSize).noalias() = (DinvBlocks[i] * res.segment(offset, blockSize)) * m_cfg.pj_relaxation;
            if (offset >= nSprings + nDampers) dlambda.segment(offset, blockSize) *= alpha;

            lambda.segment(offset, blockSize).noalias() += dlambda.segment(offset, blockSize);

            if (offset >= nSprings + nDampers) project(lambda, dlambda, mus[offset - nSprings - nDampers], offset, blockSize);

            u_corr.noalias() += MiW_T.middleCols(offset, blockSize) * dlambda.segment(offset, blockSize);

            // res_norm_sq += dlambda.segment(offset, blockSize).squaredNorm() / (DinvBlocks[i].diagonal().squaredNorm());
        }
    };

    auto pj_sweep = [&](int iter) -> void {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSweep);

        res = rhs;
        res -= W_row * u_corr;
        res -= C_vec.cwiseProduct(lambda);
        dlambda = DinvMatrix * res * m_cfg.pj_relaxation;
        if (nSprings + nDampers < numLambda) dlambda.tail(numLambda - nSprings - nDampers) *= alpha;

        lambda += dlambda;
        project_all(lambda, dlambda, mus, nSprings, nDampers, Nnfc, Nfc);
        u_corr.noalias() += MiW_T * dlambda;
    };

    if (!m_cfg.pj_nesterov) {
        VectorXr u_k = u_corr;

        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            pgs_sweep(iter);
            m_last_iters = iter + 1;

            if (check_residual(u_k, u_corr, iter) < 1.0) break;
            u_k = u_corr;
        }
    } else {
        // Nesterov acceleration outer loop
        VectorXr lambda_k = lambda;
        VectorXr u_k = u_corr;
        VectorXr lambda_y = lambda;
        VectorXr u_y = u_corr;

        double thk = 1.0;
        real_t err_prev = std::numeric_limits<real_t>::infinity();
        int restarts = 0;
        bool momentum_disabled = false;

        for (int iter = 0; iter < m_cfg.pj_max_iterations; iter++) {
            lambda = lambda_y;
            u_corr = u_y;

            pgs_sweep(iter);
            const VectorXr lambda_k1 = lambda;
            const VectorXr u_k1 = u_corr;
            m_last_iters = iter + 2;

            const real_t err = check_residual(u_y, u_k1, iter);
            if (err <= (real_t)1) break;

            const double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
            const double betak1 = (thk - 1.0) / thk1;

            bool restart = false;
            if (!std::isfinite((double)err)) {
                restart = true;
            } else if (err_prev < std::numeric_limits<real_t>::infinity() && err > (real_t)1.05 * err_prev) {
                restart = true;
            } else if (betak1 > (double)m_cfg.pj_nesterov_beta_threshold) {
                restart = true;
            } else if (!std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0) {
                restart = true;
            } else if ((u_y - u_k1).dot(u_k1 - u_k) > 0.0) {
                restart = true;
            }

            if (restart) {
                lambda_y = lambda_k1;
                u_y = u_k1;
                lambda_k = lambda_k1;
                thk = 1.0;
                if (++restarts >= m_cfg.pj_nesterov_restart_limit) momentum_disabled = true;
            } else {
                if (momentum_disabled) {
                    lambda_y = lambda_k1;
                    u_y = u_k1;
                    thk = 1.0;
                } else {
                    lambda_y = lambda_k1 + (real_t)betak1 * (lambda_k1 - lambda_k);
                    u_y = u_k1 + (real_t)betak1 * (u_k1 - u_k);
                    thk = thk1;
                }
            }
            lambda_k = lambda_k1;
            u_k = u_k1;
            err_prev = err;
        }
    }

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) WarmstartProvider::storeImpulse(lambda.segment(nSprings + nDampers, nContacts), m_dyn, /*invertNormalSign=*/true);

    return u_free - u_corr;
}

}  // namespace cardillo::solver