#include "projected_gauss_seidel.hpp"

#include <cmath>

namespace cardillo::solver {

VectorXr ProjectedGaussSeidelSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSetup);

    const int numSprings = m_dyn.numSprings();
    const int numDampers = m_dyn.numDampers();

    VectorXr u_free = m_assembler.ufree(dt, theta);

    if (numSprings + numDampers == 0) return u_free;

    const BlockDiagonal& Dinv = m_assembler.Dinv(dt, theta);
    const VectorXr rhs = m_assembler.rhs(dt, theta);
    VectorXr u_corr = VectorXr::Zero(u_free.size());

    const auto& Wg = m_dyn.Wg().asSparse();
    const auto& Wgamma = m_dyn.Wgamma().asSparse();
    const auto& C = m_dyn.Cdiag() * (1.0 / (theta * dt * dt));
    const auto& A = m_dyn.Adiag() * (1.0 / (theta * dt));

    VectorXr res = VectorXr::Zero(numSprings + numDampers);
    VectorXr lambda = VectorXr::Zero(numSprings + numDampers);

    if (m_dyn.Lambda_g().size() != numSprings) m_dyn.setLambda_g(VectorXr::Zero(numSprings));
    if (m_dyn.Lambda_gamma().size() != numDampers) m_dyn.setLambda_gamma(VectorXr::Zero(numDampers));

    if (numSprings > 0) lambda.head(numSprings) = m_dyn.Lambda_g();
    if (numDampers > 0) lambda.tail(numDampers) = m_dyn.Lambda_gamma();

    std::cout << "Setup complete, starting iterations...\n" << std::endl;

    auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidel);

    const real_t alpha = m_cfg.pj_alpha;
    for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedGaussSeidelSweep);

        u_corr = VectorXr::Zero(u_corr.size());
        if (numSprings > 0) u_corr.noalias() += Wg.transpose() * lambda.head(numSprings);
        if (numDampers > 0) u_corr.noalias() += Wgamma.transpose() * lambda.tail(numDampers);
        u_corr = m_dyn.MinvDiag().cwiseProduct(u_corr);

        res = rhs;
        if (numSprings > 0) res.head(numSprings).noalias() -= Wg * u_corr;
        if (numSprings > 0) res.head(numSprings).noalias() -= C.cwiseProduct(lambda.head(numSprings));
        if (numDampers > 0) res.tail(numDampers).noalias() -= Wgamma * u_corr;
        if (numDampers > 0) res.tail(numDampers).noalias() -= A.cwiseProduct(lambda.tail(numDampers));
        lambda.noalias() += Dinv * res * 0.9;

        real_t res_norm = res.norm();

        if (iter % 100 == 0) std::cout << "PGS iter " << iter << ", residual norm: " << res_norm << std::endl;

        if (res_norm < m_cfg.pj_tol_abs) break;
        m_last_iters = iter + 1;
    }

    if (numSprings > 0) m_dyn.setLambda_g(lambda.head(numSprings));
    if (numDampers > 0) m_dyn.setLambda_gamma(lambda.tail(numDampers));

    return u_free - u_corr;
}

}  // namespace cardillo::solver