#include "conjugate_gradient.hpp"
#include <cmath>

namespace cardillo::solver {

VectorXr ConjugateGradientSolver::solve(real_t dt, real_t theta) {
    if (m_dyn.numContacts() > 0) throw std::runtime_error("Conjugate Gradient does not handle contact constraints");

    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConjugateGradientSetup);

    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int numSprings = m_dyn.numSprings();
    const int numDampers = m_dyn.numDampers();
    if (numSprings + numDampers == 0) return u_free;

    const auto& Dinv = m_assembler.Dinv(dt, theta);

    const VectorXr rhs = m_assembler.rhs(dt, theta);
    VectorXr u_corr = VectorXr::Zero(u_free.size());

    Eigen::SparseMatrix<real_t, Eigen::RowMajor> W_row(m_assembler.W());
    Eigen::SparseMatrix<real_t, Eigen::ColMajor> W_T = W_row.transpose();  // N x M

    const VectorXr C_vec = m_assembler.C(dt, theta);
    const VectorXr Minv = m_dyn.MinvDiag();

    VectorXr res = VectorXr::Zero(numSprings + numDampers);
    VectorXr lambda = VectorXr::Zero(numSprings + numDampers);
    VectorXr real_res = VectorXr::Zero(numSprings + numDampers);

    VectorXr z = VectorXr::Zero(numSprings + numDampers);
    VectorXr p = VectorXr::Zero(numSprings + numDampers);
    VectorXr Ap = VectorXr::Zero(numSprings + numDampers);
    VectorXr u_temp = VectorXr::Zero(u_free.size());

    // Warmstart
    if (m_dyn.Lambda_g().size() != numSprings) m_dyn.setLambda_g(VectorXr::Zero(numSprings));
    if (m_dyn.Lambda_gamma().size() != numDampers) m_dyn.setLambda_gamma(VectorXr::Zero(numDampers));

    if (numSprings > 0) lambda.head(numSprings) = m_dyn.Lambda_g();
    if (numDampers > 0) lambda.tail(numDampers) = m_dyn.Lambda_gamma();

    // ---------------------- SETUP CONJUGATE GRADIENT ----------------------

    // res = b - A * x
    u_corr = Minv.cwiseProduct(W_T * lambda);
    res = rhs;
    res.noalias() -= W_row * u_corr;
    res.noalias() -= C_vec.cwiseProduct(lambda);

    z = Dinv * res;
    p = z;

    real_t rz_new = 0.0;
    real_t rz_old = z.dot(res);
    real_t res_norm, alpha, beta;

    // ---------------------- CONJUGATE GRADIENT ITERATION ----------------------

    auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConjugateGradientSolve);

    for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConjugateGradientSweep);

        // Ap = A * p
        u_temp.noalias() = W_T * p;
        u_corr.noalias() = Minv.cwiseProduct(u_temp);
        Ap.noalias() = W_row * u_corr;
        Ap.noalias() += C_vec.cwiseProduct(p);

        alpha = rz_old / p.dot(Ap);

        lambda += alpha * p;
        res -= alpha * Ap;
        z = Dinv * res;

        rz_new = z.dot(res);
        res_norm = res.norm();

        if (res_norm < m_cfg.pj_tol_abs) break;

        if (std::isnan(res_norm) || std::isinf(res_norm)) {
            std::cerr << "[ConjugateGradient] Divergence detected! res_norm = " << res_norm << "\n";
            throw std::runtime_error("Conjugate Gradient diverged");
        }

        beta = rz_new / rz_old;
        p = z + beta * p;

        rz_old = rz_new;
        m_last_iters = iter + 1;
    }

    if (m_cfg.debug_pj) {
        std::cout << "[ConjugateGradient] Finished with res_norm = " << res_norm << " after " << m_last_iters << " iterations.\n";
    }

    u_corr = Minv.cwiseProduct(W_T * lambda);

    if (numSprings > 0) m_dyn.setLambda_g(lambda.head(numSprings));
    if (numDampers > 0) m_dyn.setLambda_gamma(lambda.tail(numDampers));

    return u_free - u_corr;
}

}  // namespace cardillo::solver