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

    // Debug: dump spring/contact-level information to help diagnose no-motion
    {
        const auto &mgr = m_sys.forceManager();
        const auto &cons = mgr.constraints();
        std::cerr << "[Moreau] after refresh: numSprings=" << cons.size() << std::endl;
        for (size_t i = 0; i < cons.size(); ++i) {
            const auto &c = cons[i];
            std::cerr << "  [Spring " << i << "] xA=" << c.xA.transpose()
                      << " xB=" << c.xB.transpose()
                      << " g=" << c.g.transpose()
                      << " gdot=" << c.gdot.transpose()
                      << " K.norm=" << c.K.norm()
                      << " D.norm=" << c.D.norm()
                      << std::endl;
        }
    }

    // 4) Build rhs = h*(f_ext - dt^2 * Wg * Lambda_g)
    const auto& Wg = m_dyn.WgSparse();
    const auto& Kg = m_dyn.Kmat();
    const int Cg = (int)Wg.rows();

    // rhs = dt * f_ext - dt^2 * Wg^T * Lambda_g (Wg is Cg x totalV)
    VectorXr rhs = dt * fn;
    if (Cg > 0) {
        m_Lambda_g = (1.0 / dt) * Kg * m_dyn.gcat();
        rhs -= dt * dt * Wg.transpose() * m_Lambda_g;
    }

    std::cout << "[Moreau] Building rhs: dt=" << dt << " fn=" << fn.transpose()
              << " dt^2 * Wg^T * Lambda_g=" << (dt*dt * Wg.transpose() * m_Lambda_g).transpose()
              << " rhs=" << rhs.transpose() << std::endl;

    // 5) Solve normal impulses p using W(q_mid), G, and contact free velocity v_free
    // Debug: print sizes and vectors to help diagnose why bodies might not move
    {
        const auto& W = m_dyn.W();
        const auto& Wg2 = m_dyn.WgSparse();
        std::cerr << "[Moreau] sizes: W(rows,cols) = (" << W.rows() << "," << W.cols() << ")"
                  << " Wg(rows,cols) = (" << Wg2.rows() << "," << Wg2.cols() << ")"
                  << " fn.size=" << fn.size() << " rhs.size=" << rhs.size()
                  << " m_Lambda_g.size=" << m_Lambda_g.size() << std::endl;

        // Print a few useful vectors
        std::cerr << "  fn = " << fn.transpose() << std::endl;
        std::cerr << "  m_Lambda_g = (size=" << m_Lambda_g.size() << ") ";
        if (m_Lambda_g.size() > 0) std::cerr << m_Lambda_g.transpose();
        std::cerr << std::endl;
        std::cerr << "  rhs = " << rhs.transpose() << std::endl;
    }
    VectorXr vnp1 = m_pj.iterateWithPreliminaryVelocity(vn, rhs, m_sys.config().pj_tol_abs);

    // 6) Final position: q_{n+1} = q_mid + (h/2) * B(q_mid) * u_{n+1}
    VectorXr qnp1 = q_mid;
    integrate_quaternions(qnp1, q_mid, vnp1, offQ, offV, Nb, dt, 0.5);

    // 7) Write back final state to ECS
    m_dyn.writeStateToSystem(qnp1, vnp1);

    // 8) For no-damping tests, keep Lambda_g as diagnostic: K * g at current (already set above)
    (void)Kg; // suppress unused warning if logs disabled

    // Debug: print updated Lambda_g
    {
        std::cerr << "[Moreau] updated m_Lambda_g (size=" << m_Lambda_g.size() << ") ";
        if (m_Lambda_g.size() > 0) std::cerr << m_Lambda_g.transpose();
        std::cerr << std::endl;
    }
}

}