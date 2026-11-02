#include "naive_schur_solver.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <stdexcept>

#ifdef CARDILLO_HAVE_OPENMP
#include <omp.h>
#endif

#include "../../misc/contact_projection.hpp"
#include "../../misc/contact_rho.hpp"
#include "../../misc/graph_coloring.hpp"
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

enum class SweepMode { GaussSeidel, Jacobi, Colored };

struct SchurProblem {
    RowMatrixXXr G;  // n_c x n_c Schur complement onto the contact rows
    std::vector<ContactBlock> blocks;
    misc::Coloring coloring;
    SweepMode mode{SweepMode::GaussSeidel};
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

void gaussSeidelSweep(const SchurProblem& p, const VectorXr& c, VectorXr& pi) {
    for (const auto& blk : p.blocks) updateBlockInPlace(p, blk, blockVelocity(p, blk, c, pi), pi);
}

void jacobiSweep(const SchurProblem& p, const VectorXr& c, VectorXr& pi) {
    const VectorXr xi = c + p.G * pi;
    const int n = (int)p.blocks.size();
#pragma omp parallel for schedule(static)
    for (int k = 0; k < n; ++k) {
        const auto& blk = p.blocks[(size_t)k];
        Vector3r xiBlk = Vector3r::Zero();
        xiBlk.head(blk.dim) = xi.segment(blk.offset, blk.dim);
        updateBlockInPlace(p, blk, xiBlk, pi);
    }
}

void coloredSweep(const SchurProblem& p, const VectorXr& c, VectorXr& pi) {
    std::vector<Vector3r> xi(p.blocks.size(), Vector3r::Zero());
#pragma omp parallel
    for (const auto& colorClass : p.coloring.colorClasses) {
        const int n = (int)colorClass.size();
#pragma omp for schedule(static)
        for (int k = 0; k < n; ++k) {
            const int block = colorClass[(size_t)k];
            xi[(size_t)k] = blockVelocity(p, p.blocks[(size_t)block], c, pi);
        }
#pragma omp for schedule(static)
        for (int k = 0; k < n; ++k) {
            const int block = colorClass[(size_t)k];
            updateBlockInPlace(p, p.blocks[(size_t)block], xi[(size_t)k], pi);
        }
    }
}

void sweep(const SchurProblem& p, const VectorXr& c, VectorXr& pi) {
    switch (p.mode) {
        case SweepMode::Jacobi:
            jacobiSweep(p, c, pi);
            break;
        case SweepMode::Colored:
            coloredSweep(p, c, pi);
            break;
        default:
            gaussSeidelSweep(p, c, pi);
    }
}

misc::Coloring colorCouplingGraph(const SchurProblem& p) {
    const int n = (int)p.blocks.size();
    std::vector<std::vector<int>> adjacency((size_t)n);
    for (int i = 0; i < n; ++i) {
        const auto& bi = p.blocks[(size_t)i];
        for (int j = i + 1; j < n; ++j) {
            const auto& bj = p.blocks[(size_t)j];
            const bool coupled = !p.G.block(bi.offset, bj.offset, bi.dim, bj.dim).isZero(0) || !p.G.block(bj.offset, bi.offset, bj.dim, bi.dim).isZero(0);
            if (coupled) {
                adjacency[(size_t)i].push_back(j);
                adjacency[(size_t)j].push_back(i);
            }
        }
    }
    return misc::colorGreedyWelshPowell(n, adjacency);
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

    const std::string& sweepMode = m_cfg.condensed_sweep_mode;
    p.mode = (sweepMode == "jacobi") ? SweepMode::Jacobi : (sweepMode == "colored") ? SweepMode::Colored : SweepMode::GaussSeidel;
    if (sweepMode == "chaotic" && m_cfg.debug_pj) std::cout << "[NaiveSchur] condensed.sweep_mode=chaotic is not supported; using gauss_seidel." << std::endl;
    if (p.mode == SweepMode::Colored) p.coloring = colorCouplingGraph(p);

    p.stepScale = m_cfg.pj_alpha * m_cfg.pj_relaxation;
    p.useNewton = (m_cfg.condensed_local_solve == "newton");
    p.newton =
        NewtonACParams{m_cfg.condensed_newton_max_iters, m_cfg.condensed_newton_tol, (m_cfg.condensed_newton_rho_strategy == "full") ? NewtonRhoStrategy::FullSpectral : NewtonRhoStrategy::Split};

    VectorXr pi = VectorXr::Zero(nC);
    if (m_cfg.pj_warmstart) WarmstartProvider::applyWarmstart(pi, m_dyn);

#ifdef CARDILLO_HAVE_OPENMP
    if (m_cfg.condensed_num_threads > 0) omp_set_num_threads(m_cfg.condensed_num_threads);
#endif

    const bool useNesterov = m_cfg.pj_nesterov;
    const bool useChebyshev = m_cfg.pj_chebyshev && !useNesterov;
    const bool useAnderson = m_cfg.pj_anderson && !useNesterov && !useChebyshev;
    if (useChebyshev && m_cfg.debug_pj) std::cout << "[NaiveSchur] pj.chebyshev is not supported; using the standard sweep." << std::endl;
    if (m_cfg.debug_pj) std::cout << "[NaiveSchur] nContactRows=" << nC << " nBlocks=" << p.blocks.size() << " colors=" << p.coloring.numColors << std::endl;
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

    if (useNesterov) {
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
    } else if (useAnderson) {
        const int m = std::max(1, m_cfg.pj_anderson_m);
        std::deque<VectorXr> zHist, fHist;
        VectorXr z = pi;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            VectorXr w = z;
            doSweep(w);
            m_last_iters = iter + 1;
            if (!w.allFinite()) fail(iter);

            const real_t err = residual(velocity(z), velocity(w));
            m_lastError = err;
            if (err <= (real_t)1) {
                z = w;
                break;
            }

            const VectorXr f = w - z;
            zHist.push_back(z);
            fHist.push_back(f);
            if ((int)fHist.size() > m + 1) {
                zHist.pop_front();
                fHist.pop_front();
            }
            const int mk = (int)fHist.size() - 1;
            if (mk < 1) {
                z = w;
                continue;
            }
            MatrixXXr dF(nC, mk), dZ(nC, mk);
            for (int j = 0; j < mk; ++j) {
                dF.col(j) = fHist[(size_t)j + 1] - fHist[(size_t)j];
                dZ.col(j) = zHist[(size_t)j + 1] - zHist[(size_t)j];
            }
            const VectorXr gamma = dF.colPivHouseholderQr().solve(f);
            const VectorXr candidate = w - (dZ + dF) * gamma;
            if (gamma.allFinite() && candidate.allFinite()) {
                z = candidate;
            } else {
                z = w;
                zHist.clear();
                fHist.clear();
            }
        }
        pi = z;
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
