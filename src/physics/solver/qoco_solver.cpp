#include "qoco_solver.hpp"
#include <dlfcn.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <qoco/qoco.h>

namespace cardillo::solver {

namespace {

template <typename T>
T loadSymbol(void* handle, const char* name) {
    dlerror();
    void* sym = dlsym(handle, name);
    const char* err = dlerror();
    if (err != nullptr || sym == nullptr) {
        throw std::runtime_error(std::string("Failed to resolve symbol '") + name + "': " + (err ? err : "unknown error"));
    }
    return reinterpret_cast<T>(sym);
}

}  // namespace

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(dyn), m_cfg(cfg) {
    m_assembler = cardillo::physics::assembly::QocoAssembler(m_dyn);
}

QocoSolver::~QocoSolver() {
    if (m_qoco_solver) {
        m_qoco_cleanup(m_qoco_solver);
        m_qoco_solver = nullptr;
    }
    if (m_qoco_lib_handle) {
        dlclose(m_qoco_lib_handle);
        m_qoco_lib_handle = nullptr;
    }
}

void QocoSolver::ensureQocoApiLoaded(bool first_init) {
    if (m_qoco_lib_handle) return;

    std::vector<std::pair<std::string, std::string>> candidates;
#ifdef QOCO_CPU_LIB_PATH
    const std::string cpu_lib = QOCO_CPU_LIB_PATH;
#else
    const std::string cpu_lib;
#endif
#ifdef QOCO_CUDA_LIB_PATH
    const std::string cuda_lib = QOCO_CUDA_LIB_PATH;
#else
    const std::string cuda_lib;
#endif

    const std::string requested = m_cfg.qoco_backend;
    if (requested == "cpu") {
        candidates.emplace_back("cpu", cpu_lib);
    } else if (requested == "cuda") {
        candidates.emplace_back("cuda", cuda_lib);
    } else {
        // auto: prefer cuda, then cpu
        candidates.emplace_back("cuda", cuda_lib);
        candidates.emplace_back("cpu", cpu_lib);
    }

    std::string attempt_log;
    for (const auto& candidate : candidates) {
        if (candidate.second.empty()) {
            attempt_log += "  - " + candidate.first + ": library path not configured\n";
            continue;
        }
        void* handle = dlopen(candidate.second.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!handle) {
            const char* err = dlerror();
            attempt_log += "  - " + candidate.first + ": " + std::string(err ? err : "dlopen failed") + "\n";
            continue;
        }

        try {
            m_set_default_settings = loadSymbol<SetDefaultSettingsFn>(handle, "set_default_settings");
            m_qoco_setup = loadSymbol<QocoSetupFn>(handle, "qoco_setup");
            m_qoco_solve = loadSymbol<QocoSolveFn>(handle, "qoco_solve");
            m_qoco_cleanup = loadSymbol<QocoCleanupFn>(handle, "qoco_cleanup");
            m_qoco_lib_handle = handle;
            m_loaded_backend = candidate.first;

            if (first_init) {
                std::cout << "  Requested backend: " << requested << "\n";
                std::cout << "  Active backend: " << m_loaded_backend << "\n";
            }
            return;
        } catch (...) {
            dlclose(handle);
            throw;
        }
    }

    throw std::runtime_error("Failed to load QOCO runtime backend. Attempts:\n" + attempt_log);
}

void QocoSolver::initQocoSolver(real_t dt, real_t theta, bool first_init) {
    // TODO: test mix of frictional and frictionless contacts

    if (first_init) ensureQocoApiLoaded(true);

    const QOCOInt n = (QOCOInt)m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers();
    const QOCOInt m = (QOCOInt)m_dyn.numContactRows();
    const QOCOInt p = (QOCOInt)m_dyn.numSprings() + m_dyn.numDampers();
    const QOCOInt l = (QOCOInt)m_dyn.numFrictionlessContacts();
    const QOCOInt nsoc = (QOCOInt)m_dyn.numFrictionalContacts();

    QOCOCscMatrix *P, *A, *G;
    QOCOFloat *c, *b, *h;
    {
        auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::QocoAssembly);
        P = m_assembler.P(dt, theta);
        c = m_assembler.c(dt, theta);
        A = (p > 0) ? m_assembler.A(dt, theta) : nullptr;
        b = (p > 0) ? m_assembler.b(dt, theta) : nullptr;
        G = (m > 0) ? m_assembler.G(dt, theta) : nullptr;
        h = (m > 0) ? m_assembler.h(dt, theta) : nullptr;
    }
    std::vector<QOCOInt> qvec((size_t)nsoc, (QOCOInt)3);

    QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
    m_set_default_settings(settings);
    settings->verbose = m_cfg.debug_pj ? 1 : 0;
    settings->abstol = m_cfg.pj_tol_abs;
    settings->reltol = m_cfg.pj_tol_rel;
    settings->abstol_inacc = m_cfg.pj_tol_abs;
    settings->reltol_inacc = m_cfg.pj_tol_rel;
    // TODO: Add different regularizations for these parts in config
    // I think only the G part should be regularized due to rank deficiency
    settings->kkt_static_reg_P = 1e-30;  // no perturbation
    settings->kkt_static_reg_A = 1e-30;  // no perturbation
    settings->kkt_static_reg_G = m_cfg.ip_kkt_static_reg;
    settings->kkt_dynamic_reg = m_cfg.ip_kkt_dynamic_reg;
    // TODO: Adapt to new iterative refinement settings
    // settings->iter_ref_iters = m_cfg.ip_iter_ref_iters;
    settings->ir_tol = m_cfg.pj_tol_abs;  // default 1e-6
    settings->max_ir_iters = 10;          // default 5
    settings->max_iters = m_cfg.pj_max_iterations;
    settings->ruiz_iters = 2;  // default 0

    if (first_init) {
        std::cout << "Initializing QOCO solver with settings:\n";
        std::cout << "  Absolute tolerance: " << settings->abstol << "\n";
        std::cout << "  Relative tolerance: " << settings->reltol << "\n";
        std::cout << "  Max iterations: " << settings->max_iters << "\n";
        std::cout << "  KKT static regularization P: " << settings->kkt_static_reg_P << "\n";
        std::cout << "  KKT static regularization A: " << settings->kkt_static_reg_A << "\n";
        std::cout << "  KKT static regularization G: " << settings->kkt_static_reg_G << "\n";
        std::cout << "  KKT dynamic regularization: " << settings->kkt_dynamic_reg << "\n";
        std::cout << "  Problem dimensions: n=" << n << ", m=" << m << ", p=" << p << ", l=" << l << ", nsoc=" << nsoc << "\n";

        if (m_cfg.moreau_implicit_gyroscopy)
            std::cerr << "Warning: QOCO solver does not support implicit gyroscopic forces; "
                         "ignoring config setting.\n";

        if (m_cfg.moreau_lambda_theta)
            std::cerr << "Warning: QOCO solver does not support lambda theta integration; ignoring "
                         "config setting.\n";
    }

    auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::QocoSetup);

    m_qoco_solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));
    QOCOInt exit = m_qoco_setup(m_qoco_solver, n, m, p, P, c, A, b, G, h, l, nsoc, nsoc ? qvec.data() : nullptr, settings);

    if (exit != 0) throw std::runtime_error("Failed to initialize QOCO solver");
    std::free(settings);
}

void QocoSolver::updateQocoSolver(real_t dt, real_t theta) {
    m_qoco_cleanup(m_qoco_solver);
    initQocoSolver(dt, theta, false);

    // Not working as sparsity pattern changes with time, need to re-setup the solver
    // qoco_update_matrix_data(m_qoco_solver, 0, m_assembler.A(dt, theta), m_assembler.G(dt,
    // theta)); qoco_update_vector_data(m_qoco_solver, m_assembler.c(dt, theta), m_assembler.b(dt,
    // theta), m_assembler.h(dt, theta));
}

VectorXr QocoSolver::solve(real_t dt, real_t theta) {
    if (!m_qoco_solver)
        initQocoSolver(dt, theta);
    else
        updateQocoSolver(dt, theta);

    auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::QocoSolve);
    QOCOInt exit = m_qoco_solve(m_qoco_solver);

    if (exit != QOCO_SOLVED) {
        auto exit_str = (exit == QOCO_UNSOLVED)            ? "Unsolved"
                        : (exit == QOCO_SOLVED_INACCURATE) ? "Solved Inaccurately"
                        : (exit == QOCO_NUMERICAL_ERROR)   ? "Numerical Error"
                        : (exit == QOCO_MAX_ITER)          ? "Maximum Iterations Reached"
                                                           : "Unknown Error";
        std::cerr << "\nError solving QOCO problem: " << exit_str << " (" << exit << ")" << std::endl;

        // Print diagnostics if available
        if (m_qoco_solver && m_qoco_solver->sol) {
            QOCOSolution* sol = m_qoco_solver->sol;
            std::cerr << "  Iterations: " << sol->iters << std::endl;
            std::cerr << "  Objective value: " << sol->obj << std::endl;
            std::cerr << "  Primal residual: " << sol->pres << std::endl;
            std::cerr << "  Dual residual: " << sol->dres << std::endl;
            std::cerr << "  Duality gap: " << sol->gap << std::endl;
        }

        if (exit != QOCO_SOLVED_INACCURATE) throw std::runtime_error("QOCO solver failed to solve the problem");
    }

    QOCOSolution* sol = m_qoco_solver->sol;
    m_last_iters = sol->iters;

    // Contact Impulse
    QOCOFloat z = sol->z[0];
    VectorXr Smu = m_assembler.computeSmu();
    VectorXr impulse = VectorXr::Zero(m_dyn.numContactRows());
    for (int i = 0; i < impulse.size(); ++i) impulse[i] = sol->z[i] * Smu[i];
    cardillo::solver::WarmstartProvider::storeImpulse(impulse, m_dyn);

    // Solution vector x contains stacked [v, lambda_g, lambda_gamma]
    QOCOFloat x = sol->x[0];
    VectorXr x_vec = VectorXr::Zero(m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers());
    for (int i = 0; i < x_vec.size(); ++i) x_vec[i] = sol->x[i];

    // Track spring and damper forces, return velocity.
    const int nV = (int)m_dyn.numV();
    const int nSprings = (int)m_dyn.numSprings();
    const int nDampers = (int)m_dyn.numDampers();
    if (nSprings > 0) m_dyn.setLambda_g(x_vec.segment(nV, nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(x_vec.segment(nV + nSprings, nDampers));
    return x_vec.segment(0, nV);
}

}  // namespace cardillo::solver
