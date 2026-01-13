#include "moreau.hpp"
#include "projected_jacobi.hpp"
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
    const auto& fn_ext = m_dyn.fVecExternal();    // gravity + applied external forces
    const auto& fn_gyro = m_dyn.fVecGyroscopic(); // gyroscopic forces from current state

    // 2) Inplace midpoint position update
    m_sys.explicitPositionUpdate((1.0 - m_theta) * dt);

    // 3) Evaluate contacts at midpoint positions
    m_dyn.refreshCollisionsAndSprings(dt, m_theta);

    // 4) Build the full extended RHS and solve the extended system S * x = b
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

    // RHS: M*vn + dt*(f_ext + f_gyro)
    // External forces are integrated over the time step; gyroscopic forces are evaluated at current state
    VectorXr rhs = VectorXr::Zero((index_t)extV);
    rhs.segment(0, totalV) = M_diag.cwiseProduct(vn) + dt * (fn_ext + fn_gyro);
    if (nSprings > 0) rhs.segment(totalV, nSprings) = -(1.0 / (m_theta * dt * dt)) * Cdiag.cwiseProduct(Lambda_g) - (1.0 - m_theta) / m_theta * Wg * vn;
    if (nDampers > 0) rhs.segment(totalV + nSprings, nDampers) = -(1.0 - m_theta) / m_theta * Wgamma * vn;

    // 5) Solve extended system
    auto xnp1 = m_pj.solve(rhs, m_sys.config().pj_tol_abs);
    // auto xnp1 = m_dyn.solveS(rhs);                             // No contacts
    VectorXr vnp1 = xnp1.segment(0, totalV);
    if (nSprings > 0) m_dyn.setLambda_g(xnp1.segment(totalV, nSprings)); else m_dyn.setLambda_g(VectorXr(0));

    // 6) Write final velocity to ECS
    m_dyn.writeVelocityToSystem(vnp1);

    // 7) Inplace final position update
    m_sys.explicitPositionUpdate(m_theta * dt);
}

}