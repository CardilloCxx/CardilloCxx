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
    const bool dbg = [](){ const char* e = std::getenv("MOREAU_DEBUG"); if (!e) return false; std::string v(e); return v=="1" || v=="true"; }();
    m_dyn.refreshState();
    const auto& qn = m_dyn.qVec();
    const auto& vn = m_dyn.vVec();
    const auto& fn = m_dyn.fVec();
    const auto& MinvDiag = m_dyn.MinvDiag();
    const auto& offQ = m_dyn.bodyPosOffsets();
    const auto& offV = m_dyn.bodyVelOffsets();
    const int Nb = (int)offV.size() - 1;

    // External-only acceleration
    VectorXr a_ext = MinvDiag.cwiseProduct(fn);

    // Helper: apply B(q) * u (maps velocity block to qdot)
    auto B_times_u = [](const VectorXr& q, const VectorXr& u) -> VectorXr {
        if (q.size() == 3 && u.size() == 3) {
            // point mass: qdot = [v]
            return u;
        } else if (q.size() == 7 && u.size() == 6) {
            // rigid body: q = [x(3), w, px, py, pz], u = [v(3), omega_body(3)]
            VectorXr qdot(7); qdot.setZero();
            // position part
            qdot[0] = u[0]; qdot[1] = u[1]; qdot[2] = u[2];
            const real_t w = q[3];
            const Eigen::Vector3d p(q[4], q[5], q[6]);
            const Eigen::Vector3d om(u[3], u[4], u[5]);
            // quaternion kinematics in body frame: [w_dot; p_dot] = 0.5 * [ -p^T; w I + [p]_x ] * omega
            const real_t wdot = (real_t)0.5 * (-(p.dot(om)));
            const Eigen::Vector3d pdot = (real_t)0.5 * (w * om + p.cross(om));
            qdot[3] = wdot; qdot[4] = (real_t)pdot.x(); qdot[5] = (real_t)pdot.y(); qdot[6] = (real_t)pdot.z();
            return qdot;
        }
        // Fallback: sizes unknown
        return VectorXr::Zero(q.size());
    };

    auto normalize_quat_in_q = [](VectorXr& q) {
        if (q.size() == 7) {
            Quaternion4r qq(q[3], q[4], q[5], q[6]);
            qq.normalize();
            q[3] = qq.w(); q[4] = qq.x(); q[5] = qq.y(); q[6] = qq.z();
        }
    };

    // 2) Midpoint position and preliminary velocity (external forces only)
    VectorXr q_mid = qn;
    VectorXr v_pre = vn + dt * a_ext;
    for (int b = 0; b < Nb; ++b) {
        int q0 = offQ[(size_t)b], qnxt = offQ[(size_t)b+1]; int nq = qnxt - q0;
        int v0 = offV[(size_t)b], vnxt = offV[(size_t)b+1]; int nv = vnxt - v0;
        if (nq > 0 && nv > 0) {
            VectorXr qdot_mid = B_times_u(qn.segment(q0, nq), vn.segment(v0, nv));
            if (qdot_mid.size() == nq) {
                q_mid.segment(q0, nq) = qn.segment(q0, nq) + (real_t)0.5 * dt * qdot_mid;
                VectorXr qb = q_mid.segment(q0, nq);
                normalize_quat_in_q(qb);
                if (qb.size() == nq) q_mid.segment(q0, nq) = qb;
            }
        }
    }

    // 3) Evaluate contacts at midpoint positions: write q_mid into ECS (v stays vn)
    m_dyn.writeStateToSystem(q_mid, vn);
    m_dyn.refreshCollisions();

    // 4) Solve normal impulses p using W(q_mid), G, and preliminary velocity v_pre
    // Default: Projected Jacobi; override with env USE_CPG=1 to use Conjugate Projected Gradient
    bool useCPG = false;
    if (const char* env = std::getenv("USE_CPG")) { useCPG = (std::string(env) == "1" || std::string(env) == "true"); }
    VectorXr vnp1;
    if (useCPG) {
        cardillo::solver::CPG cpg(m_dyn);
        vnp1 = cpg.solve(v_pre, (real_t)1e-5);
    } else {
    // Use configured absolute tolerance for PJ
    vnp1 = m_pj.iterateWithPreliminaryVelocity(v_pre, m_sys.config().pj_tol_abs);
    }

    // 5) Final position: q_{n+1} = q_mid + (h/2) * B(q_mid) * u_{n+1}
    VectorXr qnp1 = q_mid;
    for (int b = 0; b < Nb; ++b) {
        int q0 = offQ[(size_t)b], qnxt = offQ[(size_t)b+1]; int nq = qnxt - q0;
        int v0 = offV[(size_t)b], vnxt = offV[(size_t)b+1]; int nv = vnxt - v0;
        if (nq > 0 && nv > 0) {
            VectorXr qdot_np1 = B_times_u(q_mid.segment(q0, nq), vnp1.segment(v0, nv));
            if (qdot_np1.size() == nq) {
                qnp1.segment(q0, nq) = q_mid.segment(q0, nq) + (real_t)0.5 * dt * qdot_np1;
                VectorXr qb = qnp1.segment(q0, nq);
                normalize_quat_in_q(qb);
                if (qb.size() == nq) qnp1.segment(q0, nq) = qb;
            }
        }
    }

    // 6) Write back final state to ECS
    m_dyn.writeStateToSystem(qnp1, vnp1);
}

}
