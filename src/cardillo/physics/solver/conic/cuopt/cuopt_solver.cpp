#include "cuopt_solver.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

#include "../../warmstart.hpp"

namespace cardillo::solver {

namespace {

cuopt::Settings makeSettings(const config::Config& cfg) {
    cuopt::Settings s;
    s.tol_abs = cfg.ip_tol_abs;
    s.tol_rel = cfg.ip_tol_rel;
    s.time_limit = cfg.cuopt_time_limit;
    s.iteration_limit = cfg.cuopt_iteration_limit;
    if (cfg.cuopt_presolve == "off")
        s.presolve = cuopt::Presolve::Off;
    else if (cfg.cuopt_presolve == "papilo")
        s.presolve = cuopt::Presolve::Papilo;
    else if (cfg.cuopt_presolve == "pslp")
        s.presolve = cuopt::Presolve::Pslp;
    s.deterministic = cfg.cuopt_deterministic;
    s.verbose = cfg.debug_pj;
    return s;
}

}  // namespace

CuoptSolver::CuoptSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_assembler(m_dyn), m_cfg(cfg), m_solver(makeSettings(cfg)) {
    // Fails early without a GPU and moves the CUDA and cuDSS start-up out of the first step.
    std::string reason;
    if (!cuopt::deviceAvailable(&reason)) throw std::runtime_error("solver.type = cuopt, but cuOpt cannot solve on this machine: " + reason);
}

CuoptSolver::~CuoptSolver() {
    if (m_inaccurate_steps > 0) std::cout << "cuOpt: " << m_inaccurate_steps << " of " << m_steps << " steps ended " << cuopt::toString(cuopt::Status::Inaccurate) << "\n";
}

VectorXr CuoptSolver::solve(real_t dt, real_t theta) {
    const CscMatrix* P = nullptr;
    VectorXr* q = nullptr;
    const CscMatrix* A = nullptr;
    VectorXr* b = nullptr;
    VectorXr Smu;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::CuoptAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
        Smu = m_assembler.computeSmu();
        m_dims.zero = m_dyn.numSprings() + m_dyn.numDampers();
        m_dims.nonneg = m_dyn.numFrictionlessContacts();
        m_dims.soc.assign(static_cast<std::size_t>(m_dyn.numFrictionalContacts()), 3);
    }

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = m_dyn.numSprings();
    const int nDampers = m_dyn.numDampers();
    const int nC = m_dyn.numContactRows();
    if (m_dims.rows() == 0) {
        // No constraint rows: the unconstrained minimizer of 1/2 x'Px + q'x with ConeProgramAssembler's
        // diagonal P, i.e. v = M^-1 (M v_old + dt f).
        m_last_iters = 0;
        const VectorXr Pdiag = P->diagonal();
        return -q->head(nV).cwiseQuotient(Pdiag.head(nV));
    }

    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::CuoptSetup);
        m_solver.setup(*P, *q, *A, *b, m_dims);
    }
    if (m_steps == 0) {
        std::cout << "Initializing cuOpt " << cuopt::Solver::version() << " solver (GPU barrier method on the Wolfe dual) with settings:\n";
        std::cout << "  tol_abs / tol_rel: " << m_cfg.ip_tol_abs << " / " << m_cfg.ip_tol_rel << ", presolve: " << m_cfg.cuopt_presolve
                  << ", deterministic cuDSS: " << (m_cfg.cuopt_deterministic ? "on" : "off") << "\n";
        std::cout << "  Cones: zero=" << m_dims.zero << " nonneg=" << m_dims.nonneg << " soc=" << m_dims.soc.size() << "\n";
    }

    cuopt::Status status;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::CuoptSolve);
        status = m_solver.solve();
    }
    ++m_steps;
    const cuopt::Info& info = m_solver.info();
    m_last_iters = info.iterations;

    if (status != cuopt::Status::Optimal) {
        if (status == cuopt::Status::Inaccurate) ++m_inaccurate_steps;
        if (status == cuopt::Status::Failed || m_cfg.debug_pj) std::cerr << "\ncuOpt: " << info.message << " (iterations " << info.iterations << ")\n";
        if (status == cuopt::Status::Failed) {
            std::string msg = "cuOpt solver failed to solve the problem: " + info.message;
            if (!m_cfg.cuopt_dump_dir.empty()) {
                const std::string path = m_cfg.cuopt_dump_dir + "/step_" + std::to_string(m_steps - 1) + ".mps";
                msg += m_solver.writeProblem(path) ? "; problem written to " + path : "; could not write " + path;
            }
            throw std::runtime_error(msg);
        }
    }

    auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::CuoptExtract);
    m_solver.recover(*A, *b);
    const VectorXr& z = m_solver.z();
    VectorXr impulse = z.tail(nC).cwiseProduct(Smu);
    WarmstartProvider::storeImpulse(impulse, m_dyn);

    // Bilateral percussions are the negated multipliers of the equality (zero-cone) rows, which lead
    // z. The primal percussions agree for compliant rows only; for perfect constraints they do not
    // enter the problem.
    const VectorXr pi_b = -z.head(nSprings + nDampers);
    if (nSprings > 0) m_dyn.setLambda_g(pi_b.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(pi_b.tail(nDampers));

    if (m_cfg.debug_pj) {
        const cuopt::KktResiduals r = cuopt::kktResiduals(*P, *q, *A, *b, m_dims, m_solver.x(), m_solver.s(), z);
        std::cout << "cuOpt step " << m_steps << ": " << cuopt::toString(status) << ", " << info.iterations << " iterations, " << info.solve_time * 1e3 << " ms; KKT residuals primal " << r.primal
                  << ", dual " << r.dual << ", cone(z) " << r.dual_cone << ", complementarity " << r.complementarity << "\n";
    }
    return m_solver.x().head(nV);
}

}  // namespace cardillo::solver
