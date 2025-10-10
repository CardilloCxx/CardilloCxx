#include "moreau.hpp"

namespace cardillo::solver {

void midpointStep(cardillo::PhysicsSystem& sys, real_t dt)
{
    // Pack state
    VectorXr qn = sys.packQ();
    VectorXr vn = sys.packV();

    // Assemble M and forces
    Eigen::SparseMatrix<real_t> M = sys.assembleMassMatrix();

    // a_n = M^{-1} f_n
    VectorXr fn = sys.assembleForceVector();

    VectorXr an(sys.numV());
    an.setZero();
    VectorXr Mdiag = VectorXr::Zero(sys.numV());
    for (int k = 0; k < M.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(M, k); it; ++it) {
            if (it.row() == it.col()) Mdiag[it.row()] = it.value();
        }
    }
    for (index_t i = 0; i < sys.numV(); ++i) an[i] = fn[i] / Mdiag[i];

    VectorXr v_mid = vn + (real_t)0.5 * dt * an;
    VectorXr q_mid = qn + (real_t)0.5 * dt * vn;
    (void)q_mid; // not used for gravity-only, kept for future forces

    // Evaluate forces at midpoint (state-independent here)
    VectorXr f_mid = sys.assembleForceVector();

    VectorXr a_mid(sys.numV());
    for (index_t i = 0; i < sys.numV(); ++i) a_mid[i] = f_mid[i] / Mdiag[i];

    VectorXr vnp1 = vn + dt * a_mid;
    VectorXr qnp1 = qn + dt * v_mid;

    sys.unpackQ(qnp1);
    sys.unpackV(vnp1);
}

}
