#include "moreau.hpp"
#include "projected_jacobi.hpp"
#include <iostream>

namespace cardillo::solver {

void MoreauSolver::stepMidpoint(real_t dt)
{
    m_dyn.refreshState();
    const VectorXr& qn = m_dyn.q();
    const VectorXr& vn = m_dyn.v();
    const VectorXr& fn = m_dyn.f();
    const VectorXr& MinvDiag = m_dyn.MinvDiag();
    const index_t V0 = m_dyn.numV();

    // External-only acceleration
    VectorXr a_ext(m_dyn.numV());
    for (index_t i = 0; i < m_dyn.numV(); ++i) a_ext[i] = MinvDiag[i] * fn[i];

    // 2) Midpoint position and preliminary velocity (external forces only)
    VectorXr q_mid = qn + (real_t)0.5 * dt * vn;
    VectorXr v_pre = vn + dt * a_ext;

    // 3) Evaluate contacts at midpoint positions: write q_mid into ECS (v stays vn)
    m_dyn.writeStateToSystem(q_mid, vn);
    m_dyn.refreshState(); // rebuild W,G based on midpoint; forces unchanged
    const index_t V1 = m_dyn.numV();

    // 4) Solve normal impulses p using W(q_mid), G, and preliminary velocity v_pre
    ProjectedJacobiSolver pjs(m_dyn);
    pjs.setAlpha((real_t)1.0);
    VectorXr p = pjs.iterateWithPreliminaryVelocity(v_pre, std::nullopt, (real_t)1e-5);
    const auto& W = m_dyn.W();

    // 5) Apply impulse: u_{n+1} = v_pre + Minv * W^T * p
    const auto& Minv = m_dyn.Minv();
    VectorXr deltaV = Minv * (W.transpose() * p);
    VectorXr vnp1 = v_pre - deltaV;

    // 6) Final position: q_{n+1} = q_mid + (h/2) * u_{n+1}
    VectorXr qnp1 = q_mid + (real_t)0.5 * dt * vnp1;

    // 7) Write back final state
    m_dyn.writeStateToSystem(qnp1, vnp1);
}

}
