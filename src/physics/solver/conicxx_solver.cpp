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

ConicxxSolver::ConicxxSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(dyn), m_assembler(m_dyn), m_cfg(cfg) {}

ConicxxSolver::~ConicxxSolver() = default;

conicxx::Settings ConicxxSolver::makeSettings() const {
    conicxx::Settings settings;
    settings.tol_feas = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.tol_gap_abs = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.tol_gap_rel = static_cast<conicxx::Scalar>(m_cfg.pj_tol_rel);
    settings.tol_infeas = static_cast<conicxx::Scalar>(m_cfg.pj_tol_abs);
    settings.max_iter = m_cfg.pj_max_iterations;
    settings.static_reg_P = static_cast<conicxx::Scalar>(m_cfg.ip_kkt_static_reg);
    settings.static_reg_A = static_cast<conicxx::Scalar>(m_cfg.ip_kkt_static_reg);
    settings.dynamic_reg_delta = static_cast<conicxx::Scalar>(m_cfg.ip_kkt_dynamic_reg);
    settings.refine_max_iter = m_cfg.ip_iter_ref_iters;
    settings.warm_start = m_cfg.conicxx_warm_start;
    settings.verbose = m_cfg.debug_pj ? 1 : 0;
    return settings;
}

bool ConicxxSolver::coneSpecChanged(const conicxx::ConeSpec& spec) const { return !sameConeSpec(spec, m_last_cone_spec); }

VectorXr ConicxxSolver::solve(real_t dt, real_t theta) {
    const conicxx::Settings settings = makeSettings();

    const SparseMatrix<Eigen::ColMajor>* P = nullptr;
    VectorXr* q = nullptr;
    const SparseMatrix<Eigen::ColMajor>* A = nullptr;
    VectorXr* b = nullptr;
    conicxx::ConeSpec cone_spec;
    {
        auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConicxxAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
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
        auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConicxxSetup);
        if (!need_setup) {
            m_solver.setSettings(settings);
            if (!m_solver.updateData(P, q, A, b)) need_setup = true;
        }
        if (need_setup) {
            if (!m_solver.setup(*P, *q, *A, *b, cone_spec, settings)) {
                throw std::runtime_error("ConicXX solver setup() failed (invalid problem data or cone spec)");
            }
            m_last_cone_spec = cone_spec;
            m_setup_done = true;

            if (!m_cfg.debug_pj) {
                // no-op; keep quiet unless debug logging requested
            } else {
                std::cout << "Initializing ConicXX solver with settings:\n";
                std::cout << "  Absolute tolerance: " << settings.tol_gap_abs << "\n";
                std::cout << "  Relative tolerance: " << settings.tol_gap_rel << "\n";
                std::cout << "  Max iterations: " << settings.max_iter << "\n";
                std::cout << "  Warm start: " << (settings.warm_start ? "on" : "off") << "\n";
                std::cout << "  Cone dims: zero=" << cone_spec.zero_dim << " nonneg=" << cone_spec.nonneg_dim << " soc_blocks=" << cone_spec.soc_dims.size() << "\n";
            }
        }
    }

    const conicxx::Solution* solution = nullptr;
    {
        auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ConicxxSolve);
        solution = &m_solver.solve();
    }
    m_last_iters = solution->info.iterations;

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
    cardillo::solver::WarmstartProvider::storeImpulse(impulse, m_dyn);

    VectorXr x_vec = VectorXr::Zero(m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers());
    const int x_copy_len = std::min<int>(x_vec.size(), static_cast<int>(solution->x.size()));
    for (int i = 0; i < x_copy_len; ++i) {
        x_vec[i] = static_cast<real_t>(solution->x[i]);
    }

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = static_cast<int>(m_dyn.numSprings());
    const int nDampers = static_cast<int>(m_dyn.numDampers());
    if (nSprings > 0) m_dyn.setLambda_g(x_vec.segment(nV, nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(x_vec.segment(nV + nSprings, nDampers));

    return x_vec.segment(0, nV);
}

}  // namespace cardillo::solver
