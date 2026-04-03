#include "qoco_solver.hpp"
#include <iostream>

#include <qoco/qoco.h>

namespace cardillo::solver {

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                       const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_cfg(cfg) {
        m_assembler = cardillo::physics::assembly::QocoAssembler(m_dyn);
    }

void QocoSolver::initQocoSolver(real_t dt, real_t theta) {
    QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
    set_default_settings(settings);
    settings->verbose = m_cfg.debug_pj ? 1 : 0;
    settings->abstol = m_cfg.pj_tol_abs;
    settings->reltol = m_cfg.pj_tol_rel;
    settings->kkt_static_reg = 1e-15;    
    settings->kkt_dynamic_reg = 1e-10;
    settings->iter_ref_iters = 10;
    settings->max_iters = m_cfg.pj_max_iterations;

    if(m_cfg.moreau_implicit_gyroscopy) 
        std::cerr << "Warning: QOCO solver does not support implicit gyroscopic forces; ignoring config setting.\n";

    if(m_cfg.moreau_lambda_theta) 
        std::cerr << "Warning: QOCO solver does not support lambda theta integration; ignoring config setting.\n";

    m_qoco_solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

    QOCOCscMatrix P = m_assembler.P(dt, theta);
    QOCOFloat* c = m_assembler.c(dt, theta);
    QOCOCscMatrix A = m_assembler.A(dt, theta);
    QOCOFloat* b = m_assembler.b(dt, theta);
    QOCOCscMatrix G = m_assembler.G(dt, theta);
    QOCOFloat* h = m_assembler.h(dt, theta);

    // second-order cones: assume friction cones of size 3
    const QOCOInt nsoc = static_cast<QOCOInt>(m_dyn.W().nRows() / 3);
    std::vector<QOCOInt> qvec((size_t)nsoc, (QOCOInt)3);

    QOCOInt exit = qoco_setup(
        m_qoco_solver,
        (QOCOInt) m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers(),  // n
        (QOCOInt) m_dyn.W().nRows(),                                       // m
        (QOCOInt) m_dyn.numSprings() + m_dyn.numDampers(),                 // p
        &P,
        c,
        &A,
        b,
        &G,
        h,
        (QOCOInt) 0, // dimension of non-negative orthant                 // l
        nsoc,
        nsoc ? qvec.data() : nullptr,                                     // q
        settings
    );

    if (exit != 0) {
        throw std::runtime_error("Failed to initialize QOCO solver");
    }
}

void QocoSolver::updateQocoSolver(real_t dt, real_t theta) {
    qoco_cleanup(m_qoco_solver);
    initQocoSolver(dt, theta);

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

        throw std::runtime_error("QOCO solver failed to solve the problem");
    }

    QOCOSolution *sol = m_qoco_solver->sol;
    m_last_iters = sol->iters;

    QOCOFloat x = sol->x[0]; 
    VectorXr x_vec = VectorXr::Zero(m_dyn.numV() + m_dyn.numSprings() + m_dyn.numDampers());
    for (int i = 0; i < x_vec.size(); ++i) x_vec[i] = sol->x[i];

    // cardillo::solver::WarmstartProvider::storeImpulse(p, m_dyn);

	// Track spring and damper forces, return velocity.
    const int nV = (int)m_dyn.numV();
	const int nSprings = (int)m_dyn.numSprings();
	const int nDampers = (int)m_dyn.numDampers();
	if (nSprings > 0) m_dyn.setLambda_g(x_vec.segment(nV, nSprings)); 
    if (nDampers > 0) m_dyn.setLambda_gamma(x_vec.segment(nV + nSprings, nDampers));
	return x_vec.segment(0, nV);
}

} // namespace cardillo::solver
