#include "projected_newton.hpp"

#include <cmath>
#include <limits>

namespace cardillo::solver {

    struct workspace {
        const int Nv, Nnfc, Nfc, Nc;
        const VectorXr y_free_star;
        const VectorXr mu;

        const CsrMatrix WTstar;
        const CsrMatrix Wstar;

        VectorXr z;
        VectorXr prev_z;
        Eigen::Map<VectorXr> velocity;
        Eigen::Map<VectorXr> impulse;
        Eigen::Map<VectorXr> prev_velocity;
        Eigen::Map<VectorXr> prev_impulse;

        VectorXr tmp_impulse;
        VectorXr tmp_velocity;

        int iter = 0;
        cardillo::misc::TimingManager* timer;

        workspace(int nv, int nnfc, int nfc, int nc, VectorXr y, CsrMatrix WTstar, CsrMatrix Wstar, VectorXr m, VectorXr p, VectorXr velocity, cardillo::misc::TimingManager* timer)
            : Nv(nv),
              Nnfc(nnfc),
              Nfc(nfc),
              Nc(nc),
              y_free_star(std::move(y)),
              WTstar(std::move(WTstar)),
              Wstar(std::move(Wstar)),
              mu(std::move(m)),
              z(VectorXr::Zero(nv + nc)),
              prev_z(VectorXr::Zero(nv + nc)),
              velocity(z.data(), nv),
              impulse(z.data() + nv, nc),
              prev_velocity(prev_z.data(), nv),
              prev_impulse(prev_z.data() + nv, nc),
              tmp_impulse(VectorXr::Zero(nc)),
              tmp_velocity(VectorXr::Zero(nv)),
              timer(timer) {
            this->impulse = p;
            this->velocity = velocity;
            this->prev_velocity = this->velocity;
            this->prev_impulse = this->impulse;
        }
    };

    static inline VectorXr rdiag_sparse(const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W, const VectorXr& MinvDiag, real_t alpha) {
        const int C = (int)W.rows();
        VectorXr R = VectorXr::Zero(C);

        for (int cid = 0; cid < C; ++cid) {
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
        ws->impulse.segment(0, ws->Nnfc) = ws->impulse.segment(0, ws->Nnfc).cwiseMax((real_t)0);

        for (int i = ws->Nnfc; i < ws->Nc; i += 3) {
            ws->impulse[i] = std::max(ws->impulse[i], (real_t)0);

            const real_t t1 = ws->impulse[i + 1], t2 = ws->impulse[i + 2];
            const real_t tnorm_squared = t1 * t1 + t2 * t2;

            if (tnorm_squared <= 0) continue;

            const real_t s = std::min<real_t>((real_t)1, (ws->mu[i] * ws->impulse[i]) / std::sqrt((double)tnorm_squared));

            ws->impulse[i + 1] = s * t1;
            ws->impulse[i + 2] = s * t2;
        }
    }

    static inline void pj_sweep(std::unique_ptr<workspace>& ws) {
        auto scope = ws->timer->scope(cardillo::misc::TimingManager::TimerId::ProjectedJacobiSweep);
        ws->prev_velocity = ws->velocity;
        ws->prev_impulse = ws->impulse;

        ws->impulse.noalias() -= ws->Wstar * ws->velocity;
        ws->impulse -= ws->y_free_star;

        project(ws);

        ws->velocity.noalias() = ws->WTstar * ws->impulse;
    }

    static inline real_t residual(const VectorXr& v, const VectorXr& v_prev, cardillo::config::Config& cfg) {
        const int Nv = (int)v.size();
        const VectorXr q = v.cwiseAbs().cwiseMax(v_prev.cwiseAbs()) * cfg.pj_tol_rel + VectorXr::Constant(Nv, cfg.pj_tol_abs);
        return 1.0 / sqrt((double)Nv) * (v - v_prev).cwiseQuotient(q).norm();
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

        VectorXr y_free_star = Wstar * u_free + Rdiag.cwiseProduct(dyn.contactVVec());
        const auto& contactsAll = dyn.contacts();

        VectorXr p = VectorXr::Zero(Nc);
        if (cfg.pj_warmstart) cardillo::solver::WarmstartProvider::applyWarmstart(p, dyn);

        VectorXr mus = VectorXr::Zero(Nc);
        for (int i = 0; i < dyn.numContacts(); ++i) {
            const auto& c = contactsAll[i];
            mus[c.impulse_base_index] = c.friction_mu;
        }

        VectorXr v = WTstar * p;

        return std::make_unique<workspace>(Nv, Nnfc, Nfc, Nc, y_free_star, WTstar, Wstar, mus, p, v, dyn.timings());
    }

    static inline void standard_loop(std::unique_ptr<workspace>& ws, cardillo::config::Config& cfg) {
        for (int iter = 0; iter < cfg.pj_max_iterations; ++iter) {
            pj_sweep(ws);
            ws->iter = iter + 1;

            if (iter % 10 == 0 && residual(ws->velocity, ws->prev_velocity, cfg) <= 1) break;
        }
    }

    static inline void nesterov_loop(std::unique_ptr<workspace>& ws, cardillo::config::Config& cfg) {
        VectorXr zk = ws->z;
        VectorXr zk1 = ws->z;
        VectorXr zy = ws->z;
        double thk = 1.0;
        real_t err_prev = std::numeric_limits<real_t>::infinity();
        int restarts = 0;
        bool momentum_disabled = false;

        for (int iter = 0; iter < cfg.pj_max_iterations; ++iter) {
            ws->z = zy;
            pj_sweep(ws);
            zk1 = ws->z;
            ws->iter = iter + 1;

            const real_t err = residual(zk1.head(ws->Nv), zk.head(ws->Nv), cfg);
            if (err <= (real_t)1) break;

            const double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
            const double betak1 = (thk - 1.0) / thk1;

            bool restart = false;
            if (!std::isfinite((double)err)) {
                restart = true;
            } else if (err_prev < std::numeric_limits<real_t>::infinity() && err > (real_t)1.05 * err_prev) {
                restart = true;
            } else if (betak1 > (double)cfg.pj_nesterov_beta_threshold) {
                restart = true;
            } else if (!std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0) {
                restart = true;
            } else {
                const VectorXr vy = zy.head(ws->Nv);
                const VectorXr vk1 = zk1.head(ws->Nv);
                const VectorXr vk = zk.head(ws->Nv);
                const double d = (vy - vk1).dot(vk1 - vk);
                if (d > 0.0) restart = true;
            }

            if (restart) {
                zy = zk1;
                zk = zk1;
                thk = 1.0;
                if (++restarts >= cfg.pj_nesterov_restart_limit) momentum_disabled = true;
            } else {
                if (momentum_disabled) {
                    zy = zk1;
                    thk = 1.0;
                } else {
                    zy = zk1 + (real_t)betak1 * (zk1 - zk);
                    thk = thk1;
                }
                zk = zk1;
            }
            err_prev = err;
        }

        ws->z = zk1;
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

            if (m_dyn.numContacts() == 0) return u_free;

            ws = build_workspace(m_dyn, m_cfg, u_free);
        }

        auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ProjectedNewton);

        if (m_cfg.pj_nesterov)
            nesterov_loop(ws, m_cfg);
        else
            standard_loop(ws, m_cfg);

        m_last_iters = ws->iter;

        cardillo::solver::WarmstartProvider::storeImpulse(ws->impulse, m_dyn);

        return u_free + ws->velocity;
    }
    }  // namespace cardillo::solver