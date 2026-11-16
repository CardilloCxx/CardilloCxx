#include "projected_gauss_seidel.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "../../misc/contact_projection.hpp"
#include "nesterov.hpp"

namespace cardillo::solver {

// Updates the block at `offset` from step = rho * xi (already relaxed/scaled) and writes the
// applied change into `dlambda`: bilateral blocks take the linear step lambda - step, contact blocks
// (n = 1: frictionless, n = 3: frictional) the negative variant lambda = -proj_C(step - lambda)
// (misc/contact_projection.hpp).
static inline void updateBlock(VectorXr& lambda, VectorXr& dlambda, bool isContact, const real_t mu, const int offset, const int n) {
    if (!isContact) {
        dlambda.segment(offset, n) = -dlambda.segment(offset, n);
        lambda.segment(offset, n) += dlambda.segment(offset, n);
        return;
    }
    if (n == 1) {
        const real_t lamNew = misc::proxNormal(dlambda[offset] - lambda[offset]);
        dlambda[offset] = lamNew - lambda[offset];
        lambda[offset] = lamNew;
        return;
    }
    Vector3r lamNew;
    misc::proxCoulomb(dlambda.segment<3>(offset) - lambda.segment<3>(offset), mu, lamNew);
    dlambda.segment<3>(offset) = lamNew - lambda.segment<3>(offset);
    lambda.segment<3>(offset) = lamNew;
}

VectorXr ProjectedGaussSeidelSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(misc::TimingManager::TimerId::ProjectedGaussSeidelSetup);

    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int nSprings = m_dyn.numSprings();
    const int nDampers = m_dyn.numDampers();
    const int nContacts = m_dyn.numContactRows();
    const int nV = u_free.size();
    const int numLambda = nSprings + nDampers + nContacts;

    if (m_cfg.solver_verbosity >= 2) {
        std::cout << "PGS solve: nSprings = " << nSprings << ", nDampers = " << nDampers << ", nContacts = " << nContacts << ", numLambda = " << numLambda << std::endl;
    }

    if (numLambda == 0) return u_free;

    const VectorXr mus = m_dyn.muVec();
    const auto& DinvMatrix = m_assembler.DinvDiag(dt, theta);
    const VectorXr rhs = m_assembler.rhs(dt, theta, u_free);
    VectorXr du = VectorXr::Zero(nV);

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

    if (m_cfg.pj_warmstart && nContacts > 0) lambda.segment(nSprings + nDampers, nContacts) = m_dyn.contactPercussions();

    auto sc_solve = m_dyn.timings()->scope(misc::TimingManager::TimerId::ProjectedGaussSeidel);

    const real_t alpha = m_cfg.pj_alpha;

    // Initial velocity correction du = Minv * W * lambda from the warmstarted lambda, u = u_free + du
    du = MiW_T * lambda;

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

    auto residual = [&](const VectorXr& du_k, const VectorXr& du_k1) -> real_t {
        const int Nv = (int)du_k.size();
        res_diff.noalias() = du_k1 - du_k;
        res_v0.noalias() = u_free + du_k;
        res_v1.noalias() = u_free + du_k1;
        res_q.array() = res_v0.cwiseAbs().cwiseMax(res_v1.cwiseAbs()).array() * m_cfg.pj_tol_rel + m_cfg.pj_tol_abs;
        return (real_t)(1.0 / std::sqrt((double)Nv) * res_diff.cwiseQuotient(res_q).norm());
    };

    // solver.criterion = natural_map: the increment of each block update divided by its step is the
    // natural map of the row at the state of that update (misc::naturalMapError2); for a spring or
    // damper row it is the row residual xi = W^T u + b + c lambda, scaled with W^T u + b. The steps
    // relaxation * alpha * Dinv of the contact rows are inverted once (0: no dynamic body).
    const bool naturalMap = m_cfg.solver_criterion == config::TerminationCriterion::NaturalMap;
    const real_t tolAbs = m_cfg.pj_tol_abs, tolRel = m_cfg.pj_tol_rel;
    VectorXr invStep = VectorXr::Zero(naturalMap ? numLambda : 0);
    if (naturalMap) {
        for (int i = 0; i < nBlocks; ++i) {
            if (offsets[i] < nSprings + nDampers) continue;
            for (int j = 0; j < DinvBlocks[i].rows(); ++j) {
                const real_t step = m_cfg.pj_relaxation * alpha * DinvBlocks[i](j, j);
                invStep[offsets[i] + j] = (step > (real_t)0 && std::isfinite(step)) ? (real_t)1 / step : (real_t)0;
            }
        }
    }
    real_t sweepError2 = 0;

    // Single sweep callable used by both standard and Nesterov loops
    auto pgs_sweep = [&]() -> void {
        auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::ProjectedGaussSeidelSweep);
        sweepError2 = 0;
        for (int i = 0; i < nBlocks; ++i) {
            const int blockSize = DinvBlocks[i].rows();
            const int offset = offsets[i];

            // Row velocity xi = xi^free + W^T du + C lambda of this block
            res.segment(offset, blockSize).noalias() = rhs.segment(offset, blockSize);
            res.segment(offset, blockSize).noalias() += W_row.middleRows(offset, blockSize) * du;
            res.segment(offset, blockSize).noalias() += C_vec.segment(offset, blockSize).cwiseProduct(lambda.segment(offset, blockSize));

            const bool isContact = offset >= nSprings + nDampers;
            if (naturalMap && !isContact) {
                for (int j = offset; j < offset + blockSize; ++j) sweepError2 = std::max(sweepError2, misc::naturalMapError2(res[j], res[j] - C_vec[j] * lambda[j], tolAbs, tolRel));
            }
            dlambda.segment(offset, blockSize).noalias() = (DinvBlocks[i] * res.segment(offset, blockSize)) * m_cfg.pj_relaxation;
            if (isContact) dlambda.segment(offset, blockSize) *= alpha;
            updateBlock(lambda, dlambda, isContact, isContact ? mus[offset - nSprings - nDampers] : (real_t)0, offset, blockSize);
            if (naturalMap && isContact) {
                const real_t invN = invStep[offset];
                if (blockSize == 1) {
                    sweepError2 = std::max(sweepError2, misc::naturalMapError2(dlambda[offset] * invN, res[offset], tolAbs, tolRel));
                } else {
                    const real_t invT = invStep[offset + 1];  // both tangential rows share one step
                    sweepError2 = std::max(sweepError2, misc::naturalMapError2(dlambda[offset] * invN, dlambda[offset + 1] * invT, dlambda[offset + 2] * invT, res[offset], res[offset + 1],
                                                                               res[offset + 2], tolAbs, tolRel));
                }
            }

            du.noalias() += MiW_T.middleCols(offset, blockSize) * dlambda.segment(offset, blockSize);
        }
    };

    // Stopping measure after a sweep from the velocity correction duIn, <= 1: converged.
    auto error = [&](const VectorXr& duIn) -> real_t { return naturalMap ? std::sqrt(sweepError2) : residual(duIn, du); };

    if (!m_cfg.pj_nesterov) {
        VectorXr du_k = du;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            pgs_sweep();
            m_last_iters = iter + 1;

            const real_t err = error(du_k);
            if (!std::isfinite(err)) {
                std::cerr << "[ProjectedGaussSeidel] Divergence detected after " << iter << " iterations. Error: " << err << std::endl;
                throw std::runtime_error("Projected Gauss-Seidel diverged");
            }
            if (m_cfg.solver_verbosity >= 2 && iter % 1000 == 0) std::cout << "PGS iter " << iter << ", error " << err << std::endl;
            if (err < (real_t)1) break;
            if (!naturalMap) du_k = du;
        }
    } else {
        const NesterovSettings settings{m_cfg.pj_max_iterations, m_cfg.pj_nesterov_beta_threshold, m_cfg.pj_nesterov_restart_limit, m_cfg.solver_verbosity, "PGS"};
        m_last_iters = nesterovIterate(du, lambda, nV, settings, pgs_sweep, error);
    }

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) m_dyn.setContactPercussions(lambda.segment(nSprings + nDampers, nContacts));

    return u_free + du;
}

}  // namespace cardillo::solver