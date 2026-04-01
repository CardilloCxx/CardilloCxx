#include "qoco_solver.hpp"
#include <iostream>

#include <qoco/qoco.h>

namespace cardillo::solver {

QocoSolver::QocoSolver(cardillo::physics::DynamicsAssembler& dyn,
                       const cardillo::config::Config& cfg)
    : m_dyn(dyn), m_cfg(cfg) {}

VectorXr QocoSolver::solve(VectorXr& rhs, real_t tol) {

    std::cout << "[QocoSolver] Solving with QOCO...\n";

    // Allocate settings.
    QOCOSettings* settings = (QOCOSettings*)malloc(sizeof(QOCOSettings));

    // Set default settings.
    set_default_settings(settings);
    settings->verbose = 1;

    // Allocate solver.
    QOCOSolver* solver = (QOCOSolver*)malloc(sizeof(QOCOSolver));

    VectorXr x = m_dyn.solveS(rhs);
    return x;
}

} // namespace cardillo::solver
