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
    settings->verbose = 0;
    settings->abstol = 1e-12;
    settings->reltol = 1e-12;
    settings->kkt_static_reg = 1e-15;
    settings->kkt_dynamic_reg = 1e-15;
    settings->iter_ref_iters = 0;

    //TODO: set settings from config        

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
    if (m_dyn.timings()) { auto sc = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::QocoSolve); (void)sc; }

    // TODO: dont build and factor S for QOCO solver
    // TODO: check that implicit gyrosopic forces are disabled for QOCO solver
      m_dyn.updateStateDependentTerms();
    if(!m_qoco_solver) initQocoSolver(dt, theta); 
    else updateQocoSolver(dt, theta);

    QOCOInt exit = qoco_solve(m_qoco_solver);

    if (exit != QOCO_SOLVED) {
        std::cerr << "Error solving QOCO problem: " << exit << std::endl;
        throw std::runtime_error("QOCO solver failed");
    }

    QOCOSolution *sol = m_qoco_solver->sol;

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
