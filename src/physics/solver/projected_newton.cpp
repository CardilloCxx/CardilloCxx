#include "projected_newton.hpp"

namespace cardillo::solver {

    struct workspace {
        const int Nv, Nnfc, Nfc, Nc;
        const VectorXr y_free_star;

        const CsrMatrix WTstar;
        const CsrMatrix Wstar;

        const VectorXr mu;
        const VectorXr zeroNnfc;

        VectorXr velocity;
        VectorXr prev_impulse;
        VectorXr impulse;

        int iter = 0;

        workspace(int nv, int nnfc, int nfc, int nc, 
                VectorXr y, CsrMatrix WTstar, CsrMatrix Wstar, VectorXr m, VectorXr z, VectorXr p, VectorXr prev_p, VectorXr velocity)
            : Nv(nv), Nnfc(nnfc), Nfc(nfc), Nc(nc), 
            y_free_star(std::move(y)), WTstar(std::move(WTstar)), Wstar(std::move(Wstar)),
            mu(std::move(m)), zeroNnfc(std::move(z)), impulse(std::move(p)), prev_impulse(std::move(prev_p)), velocity(std::move(velocity)) {}
    };

    static inline VectorXr rdiag_sparse(const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W, const VectorXr& MinvDiag, real_t alpha) {
        const int C = (int)W.rows();
        VectorXr R = VectorXr::Zero(C);

        for (int cid = 0; cid < C; ++cid) {
            if (cid < 0 || cid >= C) continue;
            real_t Dii = 0;
            for (Eigen::SparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(W, cid); it; ++it) {
                const int col = it.col();
                const real_t w = it.value();
                Dii += w * w * MinvDiag[col];
            }
            R[cid] = (Dii > (real_t)0) ? (alpha / Dii) : (real_t)0;
        }
        return R;
    }

    static inline void project(std::unique_ptr<workspace>& ws) {
            
        ws->impulse.segment(0, ws->Nnfc) = ws->impulse.segment(0, ws->Nnfc).cwiseMax(ws->zeroNnfc);

        for (int i = ws->Nnfc; i < ws->Nc; i += 3) {
            const real_t pn = std::max<real_t>(ws->impulse[i], (real_t)0);
            const real_t t1 = ws->impulse[i + 1], t2 = ws->impulse[i + 2];
            const real_t tnorm = std::sqrt((double)(t1 * t1 + t2 * t2));
            const real_t mu = ws->mu[i];
            const real_t s = (tnorm > (real_t)0) ? std::min<real_t>((real_t)1, (mu * pn) / (real_t)tnorm) : (real_t)1;
            ws->impulse[i] = pn;
            ws->impulse[i + 1] = s * t1;
            ws->impulse[i + 2] = s * t2;
        }
    }

    static inline void pj_sweep(std::unique_ptr<workspace>& ws) {
        ws->prev_impulse = ws->impulse;
        ws-> velocity = ws->velocity + ws->WTstar * (ws->impulse - ws->prev_impulse);
        ws->impulse = ws->y_free_star + ws->Wstar * ws->velocity;
        project(ws);
    }

    static inline std::unique_ptr<workspace> build_workspace(cardillo::physics::DynamicsAssembler& dyn, cardillo::config::Config& cfg, const VectorXr& u_free) {
        const int Nv = (int)dyn.numV();
        const int Nc = (int)dyn.numContactRows();
        const int Nnfc = (int)dyn.numFrictionlessContacts();
        const int Nfc = Nc - Nnfc;
        const auto& W = dyn.W().asSparseRowMajor();
        const auto& Minv = dyn.MinvDiag();

        auto Rdiag = rdiag_sparse(W, Minv, cfg.pj_alpha);

        CsrMatrix WTstar = Minv.asDiagonal() * W.transpose();
        CsrMatrix Wstar = Rdiag.asDiagonal() * W;

        VectorXr y_free_star = Wstar * u_free;
        const auto& contactsAll = dyn.contacts();
        const auto& zeroNnfc = VectorXr::Zero(Nnfc);

        VectorXr p = VectorXr::Zero(Nc);
        VectorXr prev_p = VectorXr::Zero(Nc);
        if(cfg.pj_warmstart) cardillo::solver::WarmstartProvider::applyWarmstart(p, dyn);

        VectorXr mus = VectorXr::Zero(Nc);
        for (int i = 0; i < dyn.numContacts(); ++i) {
            const auto& c = contactsAll[i];
            mus[c.impulse_base_index] = c.friction_mu;
        }

        VectorXr v = VectorXr::Zero(Nv);

        return std::make_unique<workspace>(Nv, Nnfc, Nfc, Nc, y_free_star, WTstar, Wstar, mus, zeroNnfc, p, prev_p, v);
    }

    static inline void standard_loop(std::unique_ptr<workspace>& ws, cardillo::config::Config& cfg) {

        for (int iter = 0; iter < cfg.pj_max_iterations; ++iter) {
            pj_sweep(ws);
            ws->iter = iter + 1;

            // Compute error
            const VectorXr diff = ws->impulse - ws->prev_impulse;
            const real_t absErr = diff.norm();
            const real_t relErr = absErr / ws->impulse.norm();

            std::cout << "absErr: " << absErr << ", relErr: " << relErr << std::endl;

            if (absErr < cfg.pj_tol_abs && relErr < cfg.pj_tol_rel) {
                std::cout << "[ProjectedNewton] Converged in " << ws->iter << " iterations, final abs error = " << absErr << ", final rel error = " << relErr << std::endl;
                break;
            }
        }
    }

    VectorXr ProjectedNewtonSolver::solve(real_t dt, real_t theta) {
        
        std::unique_ptr<workspace> ws;
        VectorXr x_free;
        VectorXr u_free;
        {
            auto sc_build = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedNewtonSetup);
            m_assembler.buildAndFactorS(dt, theta);
            VectorXr rhs = m_assembler.rhs(dt, theta);

            // Calculate continuous solution (without impacts)
            x_free = m_assembler.solveS(rhs);
            u_free = x_free.segment(0, m_dyn.numV());
            m_last_iters = 0;

            const int numV = m_dyn.numV();
            m_dyn.setLambda_g(x_free.segment(numV, m_dyn.numSprings()));
            m_dyn.setLambda_gamma(x_free.segment(numV + m_dyn.numSprings(), m_dyn.numDampers()));

            if (m_dyn.numContacts() == 0) 
                return u_free;

            ws = build_workspace(m_dyn, m_cfg, u_free);
        }

        auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedNewton);

        standard_loop(ws, m_cfg);
        m_last_iters = ws->iter;

        cardillo::solver::WarmstartProvider::storeImpulse(ws->impulse, m_dyn);

        return u_free + ws->velocity;
    }
}