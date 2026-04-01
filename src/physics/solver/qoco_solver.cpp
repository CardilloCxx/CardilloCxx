#include "qoco_solver.hpp"
#include <iostream>

#include <qoco/qoco.h>

namespace cardillo::solver {

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                       const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_cfg(cfg) {
        m_assembler = cardillo::physics::assembly::QocoAssembler(m_dyn);
    }

void QocoSolver::initQocoSolver(real_t dt) {
    QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
    set_default_settings(settings);
    settings->verbose = 0;

    //TODO: set settings from config        

    m_qoco_solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

    QOCOCscMatrix P = m_assembler.P(dt);
    QOCOFloat* c = m_assembler.c(dt);
    QOCOCscMatrix A = m_assembler.A(dt);
    QOCOFloat* b = m_assembler.b(dt);
    QOCOCscMatrix G = m_assembler.G(dt);
    QOCOFloat* h = m_assembler.h(dt);

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

void QocoSolver::updateQocoSolver(real_t dt) {
    qoco_cleanup(m_qoco_solver);
    initQocoSolver(dt);

    // Not working as sparsity pattern changes with time, need to re-setup the solver
    // qoco_update_matrix_data(m_qoco_solver, 0, m_assembler.A(dt), m_assembler.G(dt));
    // qoco_update_vector_data(m_qoco_solver, m_assembler.c(dt), m_assembler.b(dt), m_assembler.h(dt));
}

VectorXr QocoSolver::solve(VectorXr& rhs, real_t tol) {
    // TODO: dont build and factor S for QOCO solver
    // TODO: check that implicit gyrosopic forces are disabled for QOCO solver
    // TODO: get current dt from pipeline instead
    if(!m_qoco_solver) initQocoSolver(m_cfg.sim_dt); 
    else updateQocoSolver(m_cfg.sim_dt);

    QOCOInt exit = qoco_solve(m_qoco_solver);

    if (exit != QOCO_SOLVED) {
        std::cerr << "Error solving QOCO problem: " << exit << std::endl;
        throw std::runtime_error("QOCO solver failed");
    }

    QOCOSolution *sol = m_qoco_solver->sol;

    QOCOInt iters = sol->iters;
    QOCOFloat setup_time = sol->setup_time_sec;
    QOCOFloat solve_time = sol->solve_time_sec;

    std::cout << "[QocoSolver] QOCO solved in " << iters << " iterations, setup time: " << setup_time
              << " sec, solve time: " << solve_time << " sec\n";

    QOCOFloat x = sol->x[0]; 
    VectorXr x_vec(rhs.size());

    for (int i = 0; i < rhs.size(); ++i) {
        x_vec[i] = static_cast<real_t>(sol->x[i]);
    }
    
    return x_vec;
}

} // namespace cardillo::solver
