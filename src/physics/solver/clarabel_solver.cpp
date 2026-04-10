#include "clarabel_solver.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace cardillo::solver {

namespace {

const char* statusToString(clarabel::SolverStatus status) {
    switch (status) {
        case clarabel::SolverStatus::Unsolved: return "Unsolved";
        case clarabel::SolverStatus::Solved: return "Solved";
        case clarabel::SolverStatus::PrimalInfeasible: return "Primal Infeasible";
        case clarabel::SolverStatus::DualInfeasible: return "Dual Infeasible";
        case clarabel::SolverStatus::AlmostSolved: return "Almost Solved";
        case clarabel::SolverStatus::AlmostPrimalInfeasible: return "Almost Primal Infeasible";
        case clarabel::SolverStatus::AlmostDualInfeasible: return "Almost Dual Infeasible";
        case clarabel::SolverStatus::MaxIterations: return "Maximum Iterations";
        case clarabel::SolverStatus::MaxTime: return "Maximum Time";
        case clarabel::SolverStatus::NumericalError: return "Numerical Error";
        case clarabel::SolverStatus::InsufficientProgress: return "Insufficient Progress";
        case clarabel::SolverStatus::CallbackTerminated: return "Callback Terminated";
        default: return "Unknown";
    }
}

} // namespace

ClarabelSolver::ClarabelSolver(cardillo::physics::DynamicsAssembler& dyn,
                               const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_assembler(m_dyn), m_cfg(cfg) {}

ClarabelSolver::~ClarabelSolver() = default;

void ClarabelSolver::initSolver(real_t dt, real_t theta, bool first_init) {
    clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();
    settings.verbose = m_cfg.debug_pj;
    settings.max_iter = static_cast<uint32_t>(m_cfg.pj_max_iterations);
    settings.tol_gap_abs = static_cast<double>(m_cfg.pj_tol_abs);
    settings.tol_gap_rel = static_cast<double>(m_cfg.pj_tol_rel);
    settings.tol_feas = static_cast<double>(m_cfg.pj_tol_abs);
    settings.tol_infeas_abs = static_cast<double>(m_cfg.pj_tol_abs);
    settings.tol_infeas_rel = static_cast<double>(m_cfg.pj_tol_rel);
    settings.static_regularization_constant = static_cast<double>(m_cfg.ip_kkt_static_reg);
    settings.dynamic_regularization_delta = static_cast<double>(m_cfg.ip_kkt_dynamic_reg);
    settings.iterative_refinement_max_iter = static_cast<uint32_t>(m_cfg.ip_iter_ref_iters);

    const SparseMatrix<Eigen::ColMajor>* P = nullptr;
    VectorXr* q = nullptr;
    const SparseMatrix<Eigen::ColMajor>* A = nullptr;
    VectorXr* b = nullptr;
    const std::vector<clarabel::SupportedConeT<double>>* cones = nullptr;
    {
        auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ClarabelAssembly);
        P = &m_assembler.P(dt, theta);
        q = &m_assembler.q(dt, theta);
        A = &m_assembler.A(dt, theta);
        b = &m_assembler.b(dt, theta);
        cones = &m_assembler.cones();
    }

    if (first_init) {
        std::cout << "Initializing Clarabel solver with settings:\n";
        std::cout << "  Absolute tolerance: " << settings.tol_gap_abs << "\n";
        std::cout << "  Relative tolerance: " << settings.tol_gap_rel << "\n";
        std::cout << "  Max iterations: " << settings.max_iter << "\n";
        std::cout << "  Cones: " << cones->size() << "\n";
    }

    auto setup_scope = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ClarabelSetup);
    m_solver = std::make_unique<clarabel::DefaultSolver<double>>(*P, *q, *A, *b, *cones, settings);
}

void ClarabelSolver::updateSolver(real_t dt, real_t theta) {
    m_solver.reset();
    initSolver(dt, theta, false);
}

VectorXr ClarabelSolver::solve(real_t dt, real_t theta) {
    m_dyn.updateStateDependentTerms();
    if (!m_solver) initSolver(dt, theta, true);
    else updateSolver(dt, theta);

    auto solve_scope = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::ClarabelSolve);
    m_solver->solve();

    const auto solution = m_solver->solution();
    const auto info = m_solver->info();
    m_last_iters = static_cast<int>(solution.iterations);

    if (!(solution.status == clarabel::SolverStatus::Solved || solution.status == clarabel::SolverStatus::AlmostSolved)) {
        std::cerr << "\nError solving Clarabel problem: " << statusToString(solution.status) << "\n";
        std::cerr << "  Iterations: " << info.iterations << "\n";
        std::cerr << "  Primal residual: " << info.res_primal << "\n";
        std::cerr << "  Dual residual: " << info.res_dual << "\n";
        std::cerr << "  Gap abs: " << info.gap_abs << "\n";

        if (solution.status != clarabel::SolverStatus::AlmostSolved) {
            throw std::runtime_error("Clarabel solver failed to solve the problem");
        }
    }

    VectorXr Smu = m_assembler.computeSmu();
    VectorXr impulse = VectorXr::Zero(m_dyn.numContactRows());
    const int contact_dual_offset = static_cast<int>(m_dyn.numSprings() + m_dyn.numDampers());
    const int z_len = static_cast<int>(solution.z.size());
    for (int i = 0; i < impulse.size(); ++i) {
        const int dual_index = contact_dual_offset + i;
        if (dual_index < z_len) {
            impulse[i] = static_cast<real_t>(solution.z[dual_index]) * Smu[i];
        }
    }
    cardillo::solver::WarmstartProvider::storeImpulse(impulse, m_dyn);

    VectorXr x_vec = VectorXr::Zero(m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers());
    const int x_copy_len = std::min<int>(x_vec.size(), static_cast<int>(solution.x.size()));
    for (int i = 0; i < x_copy_len; ++i) {
        x_vec[i] = static_cast<real_t>(solution.x[i]);
    }

    const int nV = static_cast<int>(m_dyn.numV());
    const int nSprings = static_cast<int>(m_dyn.numSprings());
    const int nDampers = static_cast<int>(m_dyn.numDampers());
    if (nSprings > 0) m_dyn.setLambda_g(x_vec.segment(nV, nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(x_vec.segment(nV + nSprings, nDampers));

    return x_vec.segment(0, nV);
}

} // namespace cardillo::solver
