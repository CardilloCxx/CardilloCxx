#include "moreau.hpp"

namespace cardillo::solver {

void midpointStep(cardillo::PhysicsSystem& sys, real_t dt)
{
    // Pack state
    VectorXr qn = sys.packQ();
    VectorXr vn = sys.packV();

    // Mass diagonal and forces
    const VectorXr& Mdiag = sys.massDiagonal();
    VectorXr fn = sys.assembleForceVector();

    VectorXr an(sys.numV());
    an.setZero();
    for (index_t i = 0; i < sys.numV(); ++i) an[i] = fn[i] / Mdiag[i];

    VectorXr v_mid = vn + (real_t)0.5 * dt * an;
    VectorXr q_mid = qn + (real_t)0.5 * dt * vn;
    (void)q_mid; // not used for gravity-only, kept for future forces

    // For state-independent forces (gravity-only here), reuse fn
    VectorXr a_mid(sys.numV());
    for (index_t i = 0; i < sys.numV(); ++i) a_mid[i] = fn[i] / Mdiag[i];

    VectorXr vnp1 = vn + dt * a_mid;
    VectorXr qnp1 = qn + dt * v_mid;

    sys.unpackQ(qnp1);
    sys.unpackV(vnp1);
}

}
