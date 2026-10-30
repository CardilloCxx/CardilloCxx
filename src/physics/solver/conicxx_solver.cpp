#include "conicxx_solver.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace cardillo::solver {

namespace {

bool sameConeSpec(const conicxx::ConeSpec& a, const conicxx::ConeSpec& b) {
    return a.zero_dim == b.zero_dim && a.nonneg_dim == b.nonneg_dim && a.soc_dims == b.soc_dims;
}

const char* statusToString(conicxx::Status status) { return conicxx::toString(status); }

}  // namespace

ConicxxSolver::ConicxxSolver(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(dyn), m_cfg(cfg), m_assembler(m_dyn, cfg.conicxx_compliant_rows) {}

ConicxxSolver::~ConicxxSolver() = default;

conicxx::Settings ConicxxSolver::makeSettings() const {
    conicxx::Settings settings;
    settings.tol_feas = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.tol_gap_abs = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.tol_gap_rel = static_cast<conicxx::Scalar>(m_cfg.pj_tol_rel);
    settings.tol_infeas_abs = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.tol_infeas_rel = static_cast<conicxx::Scalar>(m_cfg.pj_tol_rel);
    settings.max_iter = m_cfg.pj_max_iterations;
    const auto static_reg = static_cast<conicxx::Scalar>(m_cfg.ip_kkt_static_reg);
    settings.regularization.static_P = static_reg;
    settings.regularization.static_nonneg = static_reg;
    settings.regularization.static_soc = static_reg;
    settings.regularization.static_zero = static_reg;
    settings.regularization.dynamic_delta = static_cast<conicxx::Scalar>(m_cfg.ip_kkt_dynamic_reg);
    settings.regularization.dynamic_eps = static_cast<conicxx::Scalar>(m_cfg.conicxx_dynamic_reg_eps);
    settings.refine_max_iter = m_cfg.ip_iter_ref_iters;
    settings.refine_reltol = static_cast<conicxx::Scalar>(m_cfg.conicxx_refine_tol);
    settings.max_step_fraction = static_cast<conicxx::Scalar>(m_cfg.conicxx_max_step_fraction);
    settings.equilibrate = m_cfg.conicxx_equilibrate;
    settings.equilibrate_max_iter = m_cfg.conicxx_equilibrate_max_iter;
    settings.equilibrate_min_scale = static_cast<conicxx::Scalar>(m_cfg.conicxx_equilibrate_min_scale);
    settings.equilibrate_max_scale = static_cast<conicxx::Scalar>(m_cfg.conicxx_equilibrate_max_scale);
    settings.warm_start = m_cfg.conicxx_warm_start;
    settings.warm_lambda = static_cast<conicxx::Scalar>(m_cfg.conicxx_warm_lambda);
    settings.validate_inputs = m_cfg.conicxx_validate_inputs;
    settings.verbose = m_cfg.conicxx_verbose;
    return settings;
}

void ConicxxSolver::setWarmStartFromState(const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const VectorXr& compliance) {
    const int nV = m_dyn.numV();
    const int nS = m_dyn.numSprings();
    const int nD = m_dyn.numDampers();
    const int nE = nS + nD;
    VectorXr pi = VectorXr::Zero(nE);
    if (m_dyn.Lambda_g().size() == nS) pi.head(nS) = m_dyn.Lambda_g();
    if (m_dyn.Lambda_gamma().size() == nD) pi.tail(nD) = m_dyn.Lambda_gamma();

    VectorXr x(A.cols());
    x.head(nV) = m_dyn.vVec();
    if (x.size() > nV) x.tail(nE) = pi;

    VectorXr impulse = VectorXr::Zero(m_dyn.numContactRows());
    solver::WarmstartProvider::applyWarmstart(impulse, m_dyn);
    VectorXr z(A.rows());
    z << -pi, impulse.cwiseQuotient(m_assembler.computeSmu());

    VectorXr s = b - A * x;
    s.head(nE) += compliance.cwiseProduct(z.head(nE));
    m_solver.setWarmStart(x, s, z);
}

bool ConicxxSolver::coneSpecChanged(const conicxx::ConeSpec& spec) const { return !sameConeSpec(spec, m_last_cone_spec); }

VectorXr ConicxxSolver::solve(real_t dt, real_t theta) {
    const conicxx::Settings settings = makeSettings();

    const SparseMatrix<Eigen::ColMajor>* P = nullptr;
    VectorXr* q = nullptr;
    const SparseMatrix<Eigen::ColMajor>* A = nullptr;
    VectorXr* b = nullptr;
    const VectorXr* compliance = nullptr;
    conicxx::ConeSpec cone_spec;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ConicxxAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
        compliance = &m_assembler.compliance(dt, theta);
        cone_spec = m_assembler.coneSpec();
    }

    // A changed contact set can change the cone composition (e.g. mu crossing
    // zero switches a row between the nonnegative orthant and an SOC block)
    // while leaving A's row count and sparsity pattern coincidentally
    // unchanged, which updateData() alone would not detect. So the cone spec
    // is compared explicitly, in addition to relying on updateData()'s own
    // sparsity-pattern check for everything else.
    bool need_setup = !m_setup_done || coneSpecChanged(cone_spec);

    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ConicxxSetup);
        if (!need_setup) {
            m_solver.setSettings(settings);
            if (!m_solver.updateData(P, q, A, b, compliance)) need_setup = true;
        }
        if (need_setup) {
            if (!m_solver.setup(*P, *q, *A, *b, cone_spec, settings, compliance)) {
                throw std::runtime_error("ConicXX solver setup() failed (invalid problem data or cone spec)");
            }
            m_last_cone_spec = cone_spec;
            m_setup_done = true;
            if (settings.warm_start && m_cfg.conicxx_warm_start_after_setup) setWarmStartFromState(*A, *b, *compliance);
        }
    }

    const conicxx::Solution* solution = nullptr;
    {
        auto sc = m_dyn.timings()->scope(misc::TimingManager::TimerId::ConicxxSolve);
        solution = &m_solver.solve();
    }
    m_last_iters = solution->info.iterations;
    if (m_cfg.debug_pj) {
        const conicxx::Info& info = solution->info;
        std::cout << "ConicXX: " << (need_setup ? "setup" : "update") << ", zero=" << cone_spec.zero_dim << " nonneg=" << cone_spec.nonneg_dim
                  << " soc=" << cone_spec.soc_dims.size() << ", " << info.iterations << " iterations, " << info.kkt_factorizations << " factorizations, "
                  << info.kkt_dynamic_reg_pivots << " dynamic pivots, " << info.kkt_refinement_steps << " refinement steps, setup " << 1e3 * info.setup_time_s << " ms, solve " << 1e3 * info.solve_time_s << " ms, "
                  << (info.warm_start_used ? "warm" : "cold") << " start (merit warm " << info.warm_start_merit << ", cold " << info.cold_start_merit << ")\n";
    }

    if (solution->status != conicxx::Status::Solved) {
        std::cerr << "\nError solving ConicXX problem: " << statusToString(solution->status) << "\n";
        std::cerr << "  Iterations: " << solution->info.iterations << "\n";
        std::cerr << "  Primal residual: " << solution->info.primal_residual << "\n";
        std::cerr << "  Dual residual: " << solution->info.dual_residual << "\n";
        std::cerr << "  Duality gap: " << solution->info.duality_gap << "\n";
        throw std::runtime_error("ConicXX solver failed to solve the problem");
    }

    VectorXr Smu = m_assembler.computeSmu();
    VectorXr impulse = VectorXr::Zero(m_dyn.numContactRows());
    const int contact_dual_offset = static_cast<int>(m_dyn.numSprings() + m_dyn.numDampers());
    const int z_len = static_cast<int>(solution->z.size());
    for (int i = 0; i < impulse.size(); ++i) {
        const int dual_index = contact_dual_offset + i;
        if (dual_index < z_len) {
            impulse[i] = static_cast<real_t>(solution->z[dual_index]) * Smu[i];
        }
    }
    solver::WarmstartProvider::storeImpulse(impulse, m_dyn);

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = static_cast<int>(m_dyn.numSprings());
    const int nDampers = static_cast<int>(m_dyn.numDampers());
    // Bilateral percussions are the negated multipliers of the zero-cone rows, which lead the dual
    // vector (in the augmented form the primal percussions agree for compliant rows only).
    VectorXr pi_b(nSprings + nDampers);
    for (int i = 0; i < pi_b.size(); ++i) pi_b[i] = -static_cast<real_t>(solution->z[i]);
    if (nSprings > 0) m_dyn.setLambda_g(pi_b.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(pi_b.tail(nDampers));

    return solution->x.head(nV).cast<real_t>();
}

}  // namespace cardillo::solver
