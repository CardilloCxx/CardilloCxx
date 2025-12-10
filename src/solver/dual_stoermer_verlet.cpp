#include "dual_stoermer_verlet.hpp"
#include "projected_jacobi.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <iostream>

namespace cardillo::solver {

void DualStoermerVerletSolver::stepMidpoint(real_t dt)
{
    // 1) Get current state vectors
    m_dyn.refreshState();
    // const auto& qn = m_dyn.qVec();
    // const auto& offQ = m_dyn.bodyPosOffsets();
    // const auto& offV = m_dyn.bodyVelOffsets();
    // const int Nb = (int)offV.size() - 1;

    // 2) Inplace midpoint position update
    m_sys.linearImplicitPositionUpdate(0.5 * dt);

    // 3) Evaluate contacts at midpoint positions (Stormer-Verlet variant with gyro-effective mass)
    m_dyn.refreshCollisionsAndSpringsStormerVerlet(dt);

    // 4) Build the full extended RHS and solve the extended system S * x = b
    const auto& vn = m_dyn.vVec();
    const auto& fn_ext = m_dyn.fVecExternal();
    const auto& Wg = m_dyn.WgSparse();
    const auto& Wgamma = m_dyn.WgammaSparse();
    const auto& M_diag = m_dyn.MDiag();
    const int totalV = (m_dyn.bodyVelOffsets().empty() ? 0 : m_dyn.bodyVelOffsets().back());
    const int nSprings = (int)m_dyn.Cdiag().size();
    const int nDampers = (int)m_dyn.Adiag().size();
    const int extV = totalV + nSprings + nDampers;
    const auto& Cdiag = m_dyn.Cdiag();

    // Lambda vectors may be uninitialized on first step; copy locally and ensure correct sizes
    VectorXr Lambda_g = m_dyn.Lambda_g();
    if ((int)Lambda_g.size() != nSprings) Lambda_g = VectorXr::Zero(nSprings);

    VectorXr rhs = VectorXr::Zero((index_t)extV);
    rhs.segment(0, totalV) =  M_diag.cwiseProduct(vn) + dt * fn_ext;
    if (nSprings > 0) rhs.segment(totalV, nSprings) = - (8.0 / (dt * dt)) * Cdiag.cwiseProduct(Lambda_g) - Wg * vn;
    if (nDampers > 0) rhs.segment(totalV + nSprings, nDampers) = - Wgamma * vn;

    // 5) Solve extended system
    auto xnp1 = m_pj.solve(rhs, m_sys.config().pj_tol_abs);
    // auto xnp1 = m_dyn.solveS(rhs);                             // No contacts
    VectorXr vnp1 = xnp1.segment(0, totalV);
    VectorXr Lambda_g_12 = xnp1.segment(totalV, nSprings);
    if (nSprings > 0) m_dyn.setLambda_g(Lambda_g_12 - Lambda_g); else m_dyn.setLambda_g(VectorXr(0));

    // 6) Write final velocity to ECS
    m_dyn.writeVelocityToSystem(vnp1);

    // 7) Inplace final position update
    m_sys.explicitPositionUpdate(0.5 * dt);
}

}