#include "qoco_solver.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <dlfcn.h>

#include <qoco/qoco.h>

namespace cardillo::solver {

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                       const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_cfg(cfg) {
        m_assembler = cardillo::physics::assembly::QocoAssembler(m_dyn);
    }

QocoSolver::~QocoSolver() {
    if (m_qoco_solver) {
        if (m_qoco_initialized) {
            qoco_cleanup(m_qoco_solver);
            m_qoco_initialized = false;
            m_qoco_solver = nullptr;
        } else {
            std::free(m_qoco_solver);
            m_qoco_solver = nullptr;
        }
    }
}

void QocoSolver::initQocoSolver(real_t dt, real_t theta, bool first_init) {
    ++m_setup_calls;

    QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
    if (!settings) {
        throw std::bad_alloc();
    }
    set_default_settings(settings);
    settings->verbose = m_cfg.debug_pj ? 1 : 0;
    settings->abstol = m_cfg.pj_tol_abs;
    settings->reltol = m_cfg.pj_tol_rel;
    settings->kkt_static_reg = 1e-8;    
    settings->kkt_dynamic_reg = 1e-8;
    settings->iter_ref_iters = 0;
    settings->max_iters = m_cfg.pj_max_iterations;

    m_qoco_solver = (QOCOSolver*)std::malloc(sizeof(QOCOSolver));
    if (!m_qoco_solver) {
        std::free(settings);
        throw std::bad_alloc();
    }

    QOCOCscMatrix* P = m_assembler.P(dt, theta);
    QOCOFloat* c = m_assembler.c(dt, theta);
    QOCOCscMatrix* A = m_assembler.A(dt, theta);
    QOCOFloat* b = m_assembler.b(dt, theta);
    QOCOCscMatrix* G = m_assembler.G(dt, theta);
    QOCOFloat* h = m_assembler.h(dt, theta);

    // TODO: test mix of frictional and frictionless contacts
    const QOCOInt n = (QOCOInt)m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers();
    const QOCOInt m = (QOCOInt)m_dyn.numContactRows();
    const QOCOInt p = (QOCOInt)m_dyn.numSprings() + m_dyn.numDampers();
    const QOCOInt l = (QOCOInt)m_dyn.numFrictionlessContacts();
    const QOCOInt nsoc = (QOCOInt)m_dyn.numFrictionalContacts();
    std::vector<QOCOInt> qvec((size_t)nsoc, (QOCOInt)3);

    const bool dims_changed = (n != m_prev_n) || (m != m_prev_m) || (p != m_prev_p);
    if (dims_changed) {
        std::cout << "[QOCO] setup#" << m_setup_calls
                  << " dims changed: n=" << n << ", m=" << m << ", p=" << p
                  << ", l=" << l << ", nsoc=" << nsoc << "\n";
        m_prev_n = n;
        m_prev_m = m;
        m_prev_p = p;
    }

    if (first_init) {
        std::cout << "Initializing QOCO solver with settings:\n";
        std::cout << "  Absolute tolerance: " << settings->abstol << "\n";
        std::cout << "  Relative tolerance: " << settings->reltol << "\n";
        std::cout << "  Max iterations: " << settings->max_iters << "\n";
        std::cout << "  KKT static regularization: " << settings->kkt_static_reg << "\n";
        std::cout << "  KKT dynamic regularization: " << settings->kkt_dynamic_reg << "\n";
        std::cout << "  Problem dimensions: n=" << n << ", m=" << m << ", p=" << p << ", l=" << l << ", nsoc=" << nsoc << "\n";

        if(m_cfg.moreau_implicit_gyroscopy) 
            std::cerr << "Warning: QOCO solver does not support implicit gyroscopic forces; ignoring config setting.\n";

        if(m_cfg.moreau_lambda_theta) 
            std::cerr << "Warning: QOCO solver does not support lambda theta integration; ignoring config setting.\n";
    }

    QOCOInt exit = qoco_setup(m_qoco_solver, n, m, p, P, c, A, b, G, h, l, nsoc, nsoc ? qvec.data() : nullptr, settings);

    // P comes from a dense diagonal vector conversion and owns its CSC arrays.
    std::free(P->p);
    std::free(P->i);
    std::free(P->x);
    delete P;

    // A/G are zero-copy CSC views into assembler-owned sparse caches.
    delete A;
    delete G;

    std::free(c);
    std::free(b);
    std::free(h);
    std::free(settings);

    if (exit != 0) {
        std::free(m_qoco_solver);
        m_qoco_solver = nullptr;
        throw std::runtime_error("Failed to initialize QOCO solver");
    }

    m_qoco_initialized = true;
}

void QocoSolver::updateQocoSolver(real_t dt, real_t theta) {
    if (m_qoco_solver) {
        if (m_qoco_initialized) {
            qoco_cleanup(m_qoco_solver);
            m_qoco_initialized = false;
            m_qoco_solver = nullptr;
        } else {
            std::free(m_qoco_solver);
            m_qoco_solver = nullptr; 
        }
    }

    initQocoSolver(dt, theta, false);

    // Not working as sparsity pattern changes with time, need to re-setup the solver
    // qoco_update_matrix_data(m_qoco_solver, 0, m_assembler.A(dt, theta), m_assembler.G(dt, theta));
    // qoco_update_vector_data(m_qoco_solver, m_assembler.c(dt, theta), m_assembler.b(dt, theta), m_assembler.h(dt, theta));
}

VectorXr QocoSolver::solve(real_t dt, real_t theta)  {
    auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::QocoSolve); 

    m_dyn.updateStateDependentTerms();
    if(!m_qoco_solver) initQocoSolver(dt, theta); 
    else updateQocoSolver(dt, theta);

    QOCOInt exit = qoco_solve(m_qoco_solver);

    if (exit != QOCO_SOLVED) {
        auto exit_str = (exit == QOCO_UNSOLVED) ? "Unsolved" :
                        (exit == QOCO_SOLVED_INACCURATE) ? "Solved Inaccurately" :
                        (exit == QOCO_NUMERICAL_ERROR) ? "Numerical Error" :
                        (exit == QOCO_MAX_ITER) ? "Maximum Iterations Reached" : "Unknown Error";
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

        if(exit != QOCO_SOLVED_INACCURATE) throw std::runtime_error("QOCO solver failed to solve the problem");
    }

    QOCOSolution *sol = m_qoco_solver->sol;
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

} // namespace cardillo::solver
