#include "qoco_solver.hpp"
#include <iostream>

#include <qoco/qoco.h>

namespace cardillo::solver {

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                       const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_cfg(cfg) {
        m_assembler = cardillo::physics::assembly::QocoAssembler(m_dyn);
    }

void QocoSolver::initQocoSolver() {
    if(!m_qoco_solver) {
        QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));
        set_default_settings(settings);
        settings->verbose = 1;

        //TODO: set settings from config        

        m_qoco_solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

        QOCOCscMatrix P = m_assembler.P();
        QOCOFloat* c = m_assembler.c();
        QOCOCscMatrix A = m_assembler.A();
        QOCOFloat* b = m_assembler.b();
        QOCOCscMatrix G = m_assembler.G();
        QOCOFloat* h = m_assembler.h();

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
}

VectorXr QocoSolver::solve(VectorXr& rhs, real_t tol) {

    std::cout << "[QocoSolver] Solving with QOCO...\n";
    initQocoSolver();

    //TODO: update QocoSolver with current problem data from assembler

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
