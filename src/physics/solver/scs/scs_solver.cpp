#include "scs_solver.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "../warmstart.hpp"

namespace cardillo::solver {

ScsSolver::ScsSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_assembler(m_dyn), m_cfg(cfg) {
    if (!m_cfg.scs_stats_csv.empty()) {
        m_stats.open(m_cfg.scs_stats_csv);
        if (!m_stats) throw std::runtime_error("Cannot open scs.stats_csv file '" + m_cfg.scs_stats_csv + "'");
        m_stats << "step,n,m,n_soc,status,iters,setup_ms,solve_ms,lin_sys_ms,cone_ms,accel_ms,aa_accepted,aa_rejected,scale_updates,scale,res_pri,res_dual,gap\n";
    }
}

ScsSolver::~ScsSolver() {
    if (m_inaccurate_steps > 0) std::cout << "SCS: " << m_inaccurate_steps << " of " << m_steps << " steps ended as '" << scs::toString(scs::Status::SolvedInaccurate) << "'\n";
}

scs::Settings ScsSolver::makeSettings() const {
    scs::Settings s;
    s.eps_abs = m_cfg.scs_eps_abs >= 0 ? m_cfg.scs_eps_abs : m_cfg.pj_tol_abs;
    s.eps_rel = m_cfg.scs_eps_rel >= 0 ? m_cfg.scs_eps_rel : m_cfg.pj_tol_rel;
    s.eps_infeas = m_cfg.scs_eps_infeas;
    s.max_iters = m_cfg.scs_max_iters;
    s.time_limit_secs = m_cfg.scs_time_limit_secs;
    s.alpha = m_cfg.scs_alpha;
    s.scale = (m_cfg.scs_carry_scale && m_last_scale > 0) ? m_last_scale : m_cfg.scs_scale;
    s.rho_x = m_cfg.scs_rho_x;
    s.acceleration_lookback = m_cfg.scs_acceleration_lookback;
    s.acceleration_interval = m_cfg.scs_acceleration_interval;
    s.adaptive_scale = m_cfg.scs_adaptive_scale;
    s.normalize = m_cfg.scs_normalize;
    s.verbose = m_cfg.debug_pj;
    return s;
}

// x0 = [v_old; Lambda_g; Lambda_gamma]. For the compliant rows, stationarity in the lambda_g /
// lambda_gamma columns (P = C/(theta dt^2), A = -C/(theta dt^2), q = 0) gives y = lambda, so their
// duals are the stored multipliers as well. Contact duals are y = impulse / S_mu (y_N = mu lambda_N),
// projected onto K since the transported impulses need not lie exactly in the cone. s0 is the cone
// projection of the slack b - A x0 (zero on the equality rows).
void ScsSolver::fillWarmStart(const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const VectorXr& Smu) {
    const int nV = static_cast<int>(m_dyn.numV());
    const int nS = m_dyn.numSprings();
    const int nD = m_dyn.numDampers();
    const int nC = m_dyn.numContactRows();

    VectorXr& x = m_solver.x();
    VectorXr& y = m_solver.y();
    VectorXr& s = m_solver.s();
    x.setZero();
    y.setZero();

    x.head(nV) = m_dyn.vVec();
    if (nS > 0 && m_dyn.Lambda_g().size() == nS) {
        x.segment(nV, nS) = m_dyn.Lambda_g();
        y.head(nS) = m_dyn.Lambda_g();
    }
    if (nD > 0 && m_dyn.Lambda_gamma().size() == nD) {
        x.segment(nV + nS, nD) = m_dyn.Lambda_gamma();
        y.segment(nS, nD) = m_dyn.Lambda_gamma();
    }

    VectorXr impulse = VectorXr::Zero(nC);
    WarmstartProvider::applyWarmstart(impulse, m_dyn);
    y.tail(nC) = impulse.cwiseQuotient(Smu);
    scs::projectOntoCones(y, m_dims);

    s = b - A * x;
    s.head(m_dims.zero).setZero();
    scs::projectOntoCones(s, m_dims);
}

void ScsSolver::writeStats() {
    if (!m_stats) return;
    const auto& in = m_solver.info();
    m_stats << m_steps << ',' << m_solver.n() << ',' << m_solver.m() << ',' << m_dims.soc.size() << ',' << static_cast<int>(in.status) << ',' << in.iterations << ',' << in.setup_ms << ','
            << in.solve_ms << ',' << in.lin_sys_ms << ',' << in.cone_ms << ',' << in.accel_ms << ',' << in.accel_accepted << ',' << in.accel_rejected << ',' << in.scale_updates << ','
            << in.scale << ',' << in.res_pri << ',' << in.res_dual << ',' << in.gap << '\n';
}

VectorXr ScsSolver::solve(real_t dt, real_t theta) {
    const SparseMatrix<Eigen::ColMajor>* P = nullptr;
    VectorXr* q = nullptr;
    const SparseMatrix<Eigen::ColMajor>* A = nullptr;
    VectorXr* b = nullptr;
    VectorXr Smu;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ScsAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
        Smu = m_assembler.computeSmu();
        m_dims.zero = m_dyn.numSprings() + m_dyn.numDampers();
        m_dims.nonneg = m_dyn.numFrictionlessContacts();
        m_dims.soc.assign(static_cast<std::size_t>(m_dyn.numFrictionalContacts()), 3);
    }

    if (m_dims.rows() == 0) {
        // No constraint rows (SCS requires m > 0): the unconstrained minimizer of 1/2 x'Px + q'x with
        // ClarabelAssembler's diagonal P, i.e. v = M^-1 (M v_old + dt f).
        m_last_iters = 0;
        const VectorXr Pdiag = P->diagonal();
        return -q->head(m_dyn.numV()).cwiseQuotient(Pdiag.head(m_dyn.numV()));
    }

    const scs::Settings settings = makeSettings();
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ScsSetup);
        if (!m_solver.setup(*P, *q, *A, *b, m_dims, settings)) throw std::runtime_error("SCS setup (scs_init) failed: invalid problem data or cone dimensions");
        if (m_cfg.scs_warm_start) fillWarmStart(*A, *b, Smu);
    }
    if (!m_cfg.scs_dump_dir.empty() && m_steps % m_cfg.scs_dump_every == 0) {
        scs::ProblemDump d;
        d.P = *P;
        d.q = *q;
        d.A = *A;
        d.b = *b;
        d.dims = m_dims;
        d.warm = m_cfg.scs_warm_start;
        d.scale = settings.scale;
        if (d.warm) {
            d.x0 = m_solver.x();
            d.y0 = m_solver.y();
            d.s0 = m_solver.s();
        }
        const std::string path = m_cfg.scs_dump_dir + "/step_" + std::to_string(m_steps) + ".bin";
        if (!scs::writeDump(path, d)) std::cerr << "SCS: could not write problem dump '" << path << "'\n";
    }
    if (m_steps == 0) {
        std::cout << "Initializing SCS solver (direct, AMD + QDLDL) with settings:\n";
        std::cout << "  eps_abs / eps_rel: " << settings.eps_abs << " / " << settings.eps_rel << "\n";
        std::cout << "  normalize: " << (settings.normalize ? "on" : "off") << ", adaptive scale: " << (settings.adaptive_scale ? "on" : "off")
                  << ", warm start: " << (m_cfg.scs_warm_start ? "on" : "off") << ", carry scale: " << (m_cfg.scs_carry_scale ? "on" : "off") << "\n";
        std::cout << "  Cones: zero=" << m_dims.zero << " nonneg=" << m_dims.nonneg << " soc=" << m_dims.soc.size() << "\n";
    }

    scs::Status status;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ScsSolve);
        status = m_solver.solve(m_cfg.scs_warm_start);
    }
    ++m_steps;
    const auto& info = m_solver.info();
    m_last_iters = info.iterations;
    writeStats();

    if (status != scs::Status::Solved) {
        if (status == scs::Status::SolvedInaccurate) ++m_inaccurate_steps;
        if (!scs::isAcceptable(status) || m_cfg.debug_pj) {
            std::cerr << "\nSCS: " << scs::toString(status) << " (iterations " << info.iterations << ", res_pri " << info.res_pri << ", res_dual " << info.res_dual << ", gap " << info.gap
                      << ")\n";
        }
        if (!scs::isAcceptable(status)) throw std::runtime_error(std::string("SCS solver failed to solve the problem: ") + scs::toString(status));
    }
    if (info.scale > 0) m_last_scale = info.scale;

    auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ScsExtract);
    const VectorXr& x = m_solver.x();
    const int nC = m_dyn.numContactRows();
    VectorXr impulse = m_solver.y().tail(nC).cwiseProduct(Smu);
    WarmstartProvider::storeImpulse(impulse, m_dyn);

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = m_dyn.numSprings();
    const int nDampers = m_dyn.numDampers();
    if (nSprings > 0) m_dyn.setLambda_g(x.segment(nV, nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(x.segment(nV + nSprings, nDampers));

    return x.head(nV);
}

}  // namespace cardillo::solver
