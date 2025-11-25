#include "moreau.hpp"
#include "projected_jacobi.hpp"
#include "CPG.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <iostream>

namespace cardillo::solver {

void MoreauSolver::stepMidpoint(real_t dt)
{
    // 1) Get current state vectors
    m_dyn.refreshState();
    // const auto& qn = m_dyn.qVec();
    const auto& vn = m_dyn.vVec();
    const auto& fn = m_dyn.fVec();
    // const auto& offQ = m_dyn.bodyPosOffsets();
    // const auto& offV = m_dyn.bodyVelOffsets();
    // const int Nb = (int)offV.size() - 1;

    // 2) Inplace midpoint position update
    m_sys.explicitPositionUpdate(0.5 * dt);

    // 3) Evaluate contacts at midpoint positions
    m_dyn.refreshCollisionsAndSprings(dt);

    // 4) Build the full extended RHS and solve the extended system S * x = b
    const auto& Wg = m_dyn.WgSparse();
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
    rhs.segment(0, totalV) =  M_diag.cwiseProduct(vn) + dt * fn;
    if (nSprings > 0) rhs.segment(totalV, nSprings) = - (1.0 / (dt * dt)) * Cdiag.cwiseProduct(Lambda_g);

    // 5) Solve extended system
    auto xnp1 = m_pj.solve(rhs, m_sys.config().pj_tol_abs);
    // auto xnp1 = m_dyn.solveS(rhs);                             // No contacts
    VectorXr vnp1 = xnp1.segment(0, totalV);
    if (nSprings > 0) m_dyn.setLambda_g(xnp1.segment(totalV, nSprings)); else m_dyn.setLambda_g(VectorXr(0));

    // 6) Write final velocity to ECS
    m_dyn.writeVelocityToSystem(vnp1);

    // 7) Inplace final position update
    m_sys.explicitPositionUpdate(0.5 * dt);
}

}