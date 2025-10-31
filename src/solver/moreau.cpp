#include "moreau.hpp"
#include "projected_jacobi.hpp"
#include "CPG.h"
#include "../physics/force_interaction.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <iostream>

namespace cardillo::solver {

inline VectorXr B_times_u(const VectorXr& q, const VectorXr& u)
{
    const int nq = static_cast<int>(q.size());
    const int nu = static_cast<int>(u.size());
    VectorXr qdot = VectorXr::Zero(nq);

    // Position derivative (common to both)
    qdot.head(3) = u.head(3);

    // If not a rigid body, return here
    if (nq != 7 || nu != 6) return qdot;

    const real_t w = q(6);
    const Vector3r p = q.segment<3>(3);
    const Vector3r omega = u.tail<3>();

    // dq/dt = 0.5 * [ w*I + [p]_x ; -p^T ] * ω
    qdot.segment<3>(3) = 0.5 * (w * omega + p.cross(omega));
    qdot(6) = -0.5 * p.dot(omega);

    return qdot;
}

inline void integrate_quaternions(
    VectorXr& q_out,
    const VectorXr& q_in,
    const VectorXr& v_in,
    const std::vector<int>& offQ,
    const std::vector<int>& offV,
    int Nb,
    real_t dt,
    real_t scale)
{
    for (int b = 0; b < Nb; ++b) {
        const int q0 = offQ[b], qnxt = offQ[b + 1];
        const int v0 = offV[b], vnxt = offV[b + 1];
        const int nq = qnxt - q0, nv = vnxt - v0;
        if (nq == 0 || nv == 0) continue;

        auto q_seg = q_in.segment(q0, nq);
        auto v_seg = v_in.segment(v0, nv);
        auto q_out_seg = q_out.segment(q0, nq);

        const auto qdot = B_times_u(q_seg, v_seg);
        if (qdot.size() != nq) continue;

        q_out_seg = q_seg + scale * dt * qdot;
        if (nq == 7) q_out_seg.tail<4>().normalize();
    }
}

void MoreauSolver::stepMidpoint(real_t dt)
{
    // 1) Get current state vectors
    m_dyn.refreshState();
    const auto& qn = m_dyn.qVec();
    const auto& vn = m_dyn.vVec();
    const auto& fn = m_dyn.fVec();
    const auto& offQ = m_dyn.bodyPosOffsets();
    const auto& offV = m_dyn.bodyVelOffsets();
    const int Nb = (int)offV.size() - 1;

    // 2) Midpoint position and preliminary velocity (external forces only)
    VectorXr q_mid = qn;
    integrate_quaternions(q_mid, qn, vn, offQ, offV, Nb, dt, 0.5);

    // 3) Evaluate contacts at midpoint positions: write q_mid into ECS (v stays vn)
    m_dyn.writeStateToSystem(q_mid, vn);
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

    // 6) Final position: q_{n+1} = q_mid + (h/2) * B(q_mid) * u_{n+1}
    VectorXr qnp1 = q_mid;
    integrate_quaternions(qnp1, q_mid, vnp1, offQ, offV, Nb, dt, 0.5);

    // 7) Write back final state to ECS
    m_dyn.writeStateToSystem(qnp1, vnp1);
}

}