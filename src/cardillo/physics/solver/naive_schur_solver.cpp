#include "naive_schur_solver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "../../misc/contact_projection.hpp"
#include "../../misc/contact_rho.hpp"
#include "local_contact_newton.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

namespace {

using RowMatrixXXr = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// One contact (1 row frictionless, 3 rows frictional) of the Schur-complement problem.
struct ContactBlock {
    int offset{0};
    int dim{1};
    real_t mu{0};
    Matrix33r G{Matrix33r::Zero()};  // diagonal block of the Schur complement (zero-padded for dim 1)
    Vector3r rho{Vector3r::Zero()};  // step lengths (rhoN, rhoT, rhoT)
};

struct SchurProblem {
    RowMatrixXXr G;  // n_c x n_c Schur complement onto the contact rows
    std::vector<ContactBlock> blocks;
    real_t stepScale{(real_t)1};  // pj.alpha * pj.relaxation
    bool useNewton{false};
    NewtonACParams newton{};
};

misc::ContactRhoStrategy rhoStrategy(const std::string& s) {
    if (s == "split") return misc::ContactRhoStrategy::Split;
    if (s == "full") return misc::ContactRhoStrategy::FullSpectral;
    return misc::ContactRhoStrategy::IsotropicDiagonal;
}

// Contact velocity xi_i = c_i + G_i pi of one block, zero-padded to three entries.
Vector3r blockVelocity(const SchurProblem& p, const ContactBlock& blk, const VectorXr& c, const VectorXr& pi) {
    Vector3r xi = Vector3r::Zero();
    xi.head(blk.dim) = c.segment(blk.offset, blk.dim) + p.G.middleRows(blk.offset, blk.dim) * pi;
    return xi;
}

// Local update lambda = -proj_C(rho * xi - lambda) (misc/contact_projection.hpp)
void updateBlock(const SchurProblem& p, const ContactBlock& blk, const Vector3r& xi, Vector3r& lam) {
    if (p.useNewton && blk.dim == 3) {
        Vector3r lamNewton = lam;
        if (solveContactBlockNewtonAC(blk.G, xi, blk.mu, lamNewton, p.newton)) {
            lam = lamNewton;
            return;
        }
    }
    const Vector3r step = p.stepScale * blk.rho.cwiseProduct(xi);
    if (blk.dim == 1)
        lam[0] = misc::proxNormal(step[0] - lam[0]);
    else
        misc::proxCoulomb(step - lam, blk.mu, lam);
}

void updateBlockInPlace(const SchurProblem& p, const ContactBlock& blk, const Vector3r& xi, VectorXr& pi) {
    Vector3r lam = Vector3r::Zero();
    lam.head(blk.dim) = pi.segment(blk.offset, blk.dim);
    updateBlock(p, blk, xi, lam);
    pi.segment(blk.offset, blk.dim) = lam.head(blk.dim);
}

void sweep(const SchurProblem& p, const VectorXr& c, VectorXr& pi) {
    for (const auto& blk : p.blocks) updateBlockInPlace(p, blk, blockVelocity(p, blk, c, pi), pi);
}

}  // namespace

VectorXr NaiveSchurSolver::solve(real_t dt, real_t theta) {
    auto scSetup = m_dyn.timings()->scope(misc::TimingManager::TimerId::NaiveSchurSetup);

    if (!m_assembler.buildAndFactorS(dt, theta, m_cfg.moreau_implicit_gyroscopy)) throw std::runtime_error("NaiveSchurSolver: factorization of the saddle-point matrix failed");

    const int nV = (int)m_dyn.numV();
    const int nS = m_dyn.numSprings();
    const int nD = m_dyn.numDampers();
    const int nC = m_dyn.numContactRows();
    const int nNfc = m_dyn.numFrictionlessContacts();

    const VectorXr xFree = m_assembler.solveS(m_assembler.rhs(dt, theta));
    m_last_iters = 0;
    m_lastError = 0;

    if (nC == 0) {
        m_dyn.setLambda_g(xFree.segment(nV, nS));
        m_dyn.setLambda_gamma(xFree.segment(nV + nS, nD));
        return xFree.head(nV);
    }

    // X = S^{-1} [W_c^T; 0; 0]: response of velocities and bilateral percussions to unit contact
    // percussions.
    const auto& Wc = m_dyn.W().asSparseRowMajor();
    MatrixXXr B = MatrixXXr::Zero(xFree.size(), nC);
    B.topRows(nV) = MatrixXXr(Wc.transpose());
    const MatrixXXr X = m_assembler.solveSColumns(B);

    SchurProblem p;
    p.G = Wc * X.topRows(nV);

    const auto& contactBias = m_dyn.contactVVec();
    const VectorXr c = Wc * xFree.head(nV) + contactBias + m_dyn.restitutionVec().cwiseProduct(Wc * m_dyn.vVec() + contactBias);

    const auto strategy = rhoStrategy(m_cfg.condensed_projection_rho_strategy);
    const auto& mu = m_dyn.muVec();
    for (int row = 0; row < nC;) {
        ContactBlock blk;
        blk.offset = row;
        blk.dim = (row < nNfc) ? 1 : 3;
        blk.mu = mu[row];
        blk.G.topLeftCorner(blk.dim, blk.dim) = p.G.block(row, row, blk.dim, blk.dim);
        if (blk.dim == 1) {
            blk.rho[0] = (blk.G(0, 0) > (real_t)0) ? (real_t)1 / blk.G(0, 0) : (real_t)0;
        } else {
            blk.rho = misc::computeContactRho(blk.G, strategy);
            if (!(blk.rho.array() > (real_t)0).all()) blk.rho = misc::computeContactRho(blk.G, misc::ContactRhoStrategy::IsotropicDiagonal);
            blk.rho = blk.rho.cwiseMax((real_t)0);
        }
        p.blocks.push_back(blk);
        row += blk.dim;
    }

    p.stepScale = m_cfg.pj_alpha * m_cfg.pj_relaxation;
    p.useNewton = (m_cfg.condensed_local_solve == "newton");
    p.newton =
        NewtonACParams{m_cfg.condensed_newton_max_iters, m_cfg.condensed_newton_tol, (m_cfg.condensed_newton_rho_strategy == "full") ? NewtonRhoStrategy::FullSpectral : NewtonRhoStrategy::Split};

    VectorXr pi = VectorXr::Zero(nC);
    if (m_cfg.pj_warmstart) WarmstartProvider::applyWarmstart(pi, m_dyn);

    if (m_cfg.debug_pj) std::cout << "[NaiveSchur] nContactRows=" << nC << " nBlocks=" << p.blocks.size() << std::endl;
    scSetup.~Scope();
    auto scSolve = m_dyn.timings()->scope(misc::TimingManager::TimerId::NaiveSchur);

    const auto Xu = X.topRows(nV);
    auto velocity = [&](const VectorXr& percussions) -> VectorXr { return xFree.head(nV) + Xu * percussions; };
    auto residual = [&](const VectorXr& uPrev, const VectorXr& uCur) -> real_t {
        const VectorXr q = (uPrev.cwiseAbs().cwiseMax(uCur.cwiseAbs()).array() * m_cfg.pj_tol_rel + m_cfg.pj_tol_abs).matrix();
        return (real_t)(1.0 / std::sqrt((double)nV) * (uCur - uPrev).cwiseQuotient(q).norm());
    };
    auto doSweep = [&](VectorXr& percussions) {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::NaiveSchurSweep);
        sweep(p, c, percussions);
    };
    auto fail = [&](int iter) {
        std::cerr << "[NaiveSchur] Divergence detected after " << iter << " iterations." << std::endl;
        throw std::runtime_error("NaiveSchur solver diverged");
    };

    if (m_cfg.pj_nesterov) {
        // Nesterov momentum with adaptive restart, as in the projected Jacobi and condensed solvers.
        VectorXr piK = pi, piY = pi, uK = velocity(pi), uY = uK;
        double thk = 1.0;
        real_t errPrev = std::numeric_limits<real_t>::infinity();
        int restarts = 0;
        bool momentumDisabled = false;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            pi = piY;
            doSweep(pi);
            const VectorXr uK1 = velocity(pi);
            m_last_iters = iter + 1;
            if (!pi.allFinite() || !uK1.allFinite()) fail(iter);

            const real_t err = residual(uY, uK1);
            m_lastError = err;
            if (err <= (real_t)1) break;

            const double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
            const double betak1 = (thk - 1.0) / thk1;
            const bool restart = !std::isfinite((double)err) || (std::isfinite(errPrev) && err > (real_t)1.05 * errPrev) || betak1 > (double)m_cfg.pj_nesterov_beta_threshold ||
                                 !std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0 || (uY - uK1).dot(uK1 - uK) > 0.0;

            if (restart) {
                piY = pi;
                uY = uK1;
                thk = 1.0;
                if (++restarts >= m_cfg.pj_nesterov_restart_limit) momentumDisabled = true;
            } else if (momentumDisabled) {
                piY = pi;
                uY = uK1;
                thk = 1.0;
            } else {
                piY = pi + (real_t)betak1 * (pi - piK);
                uY = uK1 + (real_t)betak1 * (uK1 - uK);
                thk = thk1;
            }
            piK = pi;
            uK = uK1;
            errPrev = err;
        }
    } else {
        VectorXr uPrev = velocity(pi);
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            doSweep(pi);
            const VectorXr u = velocity(pi);
            m_last_iters = iter + 1;
            const real_t err = residual(uPrev, u);
            m_lastError = err;
            if (!std::isfinite(err)) fail(iter);
            if (err <= (real_t)1) break;
            uPrev = u;
        }
    }

    const VectorXr x = xFree + X * pi;
    m_dyn.setLambda_g(x.segment(nV, nS));
    m_dyn.setLambda_gamma(x.segment(nV + nS, nD));
    WarmstartProvider::storeImpulse(pi, m_dyn);
    return x.head(nV);
}

}  // namespace cardillo::solver
