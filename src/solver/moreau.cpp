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
    m_dyn.refreshState();
    const auto& qnBlocks = m_dyn.q();
    const auto& vnBlocks = m_dyn.v();
    const auto& fnBlocks = m_dyn.f();
    const auto& MinvBlocks = m_dyn.MinvBlocks();
    const int Nb = (int)m_dyn.v().size();

    // External-only acceleration
    std::vector<VectorXr> a_ext_blocks(Nb);
    for (int b = 0; b < Nb; ++b) a_ext_blocks[(size_t)b] = MinvBlocks[(size_t)b] * fnBlocks[(size_t)b];

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
    std::vector<VectorXr> q_mid_blocks(Nb), v_pre_blocks(Nb);
    for (int b = 0; b < Nb; ++b) {
        const auto& qb = qnBlocks[(size_t)b];
        const auto& ub = vnBlocks[(size_t)b];
        VectorXr qdot_mid = B_times_u(qb, ub);
        q_mid_blocks[(size_t)b] = qb + (real_t)0.5 * dt * qdot_mid;
        normalize_quat_in_q(q_mid_blocks[(size_t)b]);
        v_pre_blocks[(size_t)b] = ub + dt * a_ext_blocks[(size_t)b];
    }

    // 3) Evaluate contacts at midpoint positions: write q_mid into ECS (v stays vn)
    m_dyn.writeStateToSystem(q_mid_blocks, vnBlocks);
    m_dyn.refreshCollisions();

    // 4) Solve normal impulses p using W(q_mid), G, and preliminary velocity v_pre
    // Default: Projected Jacobi; override with env USE_CPG=1 to use Conjugate Projected Gradient
    bool useCPG = false;
    if (const char* env = std::getenv("USE_CPG")) { useCPG = (std::string(env) == "1" || std::string(env) == "true"); }
    std::vector<VectorXr> vnp1_blocks;
    if (useCPG) {
        cardillo::solver::CPG cpg(m_dyn);
        vnp1_blocks = cpg.solve(v_pre_blocks, (real_t)1e-5);
    } else {
    vnp1_blocks = m_pj.iterateWithPreliminaryVelocity(v_pre_blocks, (real_t)1e-7);
    }

    // 5) Final position: q_{n+1} = q_mid + (h/2) * B(q_mid) * u_{n+1}
    std::vector<VectorXr> qnp1_blocks(Nb);
    for (int b = 0; b < Nb; ++b) {
        VectorXr qdot_np1 = B_times_u(q_mid_blocks[(size_t)b], vnp1_blocks[(size_t)b]);
        qnp1_blocks[(size_t)b] = q_mid_blocks[(size_t)b] + (real_t)0.5 * dt * qdot_np1;
        normalize_quat_in_q(qnp1_blocks[(size_t)b]);
    }

    // Optional diagnostics: verify contact resolution and torque/linear mapping
    if (const char* dbg = std::getenv("RB_DEBUG")) {
        if (std::string(dbg) == "1" || std::string(dbg) == "true") {
            const auto& Wblocks = m_dyn.WBlocks();
            const auto& blockToBody = m_dyn.WBlockToBodyAll();
            const auto& contactToBlocks = m_dyn.WBlocksFromContactAll();
            const auto& reg = m_dyn.system().ecs();
            const auto& contacts = m_dyn.contacts();
            const auto& dynToOrig = m_dyn.dynamicContactToOriginalAll();
            // Compute W*v for v_pre and v_np1
            auto contactVel = [&](const std::vector<VectorXr>& v){
                std::vector<real_t> y; y.resize(contactToBlocks.size());
                for (size_t cid = 0; cid < contactToBlocks.size(); ++cid) {
                    auto [idxA, idxB] = contactToBlocks[cid];
                    real_t s = 0;
                    if (idxA >= 0) s += (Wblocks[(size_t)idxA].row(0) * v[(size_t)blockToBody[(size_t)idxA]])(0,0);
                    if (idxB >= 0) s += (Wblocks[(size_t)idxB].row(0) * v[(size_t)blockToBody[(size_t)idxB]])(0,0);
                    y[cid] = s;
                }
                return y;
            };
            // Compute physical normal velocity using ECS state at q_mid and the provided velocity blocks
            auto contactVelPhys = [&](const std::vector<VectorXr>& v){
                std::vector<real_t> y; y.resize(contactToBlocks.size());
                for (size_t cid = 0; cid < contactToBlocks.size(); ++cid) {
                    // Map dynamic contact id -> original contact index in contacts vector
                    const auto& c = contacts[(size_t)dynToOrig[cid]];
                    const Vector3r n = c.normal; // from a -> b

                    auto getBodyLinAng = [&](entt::entity e, Vector3r& vlin, Vector3r& omegaw, Vector3r& r){
                        vlin.setZero(); omegaw.setZero(); r.setZero();
                        if (!reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(e)) return; // static body
                        const int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(e).b;
                        if (b < 0 || b >= (int)v.size()) return;
                        const auto& xb = reg.get<cardillo::PhysicsSystem::C_Position3>(e).value;
                        r = c.point - xb; // world
                        const VectorXr& vb = v[(size_t)b];
                        if (vb.size() >= 3) {
                            vlin = Vector3r(vb[0], vb[1], vb[2]);
                        }
                        if (vb.size() >= 6 && reg.any_of<cardillo::PhysicsSystem::C_Orientation>(e)) {
                            Vector3r omega_body(vb[3], vb[4], vb[5]);
                            // rotate to world: omega_world = q * omega_body
                            const Quaternion4r q = reg.get<cardillo::PhysicsSystem::C_Orientation>(e).q;
                            omegaw = q * omega_body;
                        }
                    };

                    Vector3r vA, vB, wAw, wBw, rA, rB; vA.setZero(); vB.setZero(); wAw.setZero(); wBw.setZero(); rA.setZero(); rB.setZero();
                    getBodyLinAng(c.a, vA, wAw, rA);
                    getBodyLinAng(c.b, vB, wBw, rB);
                    // y = n · (vB + wB×rB − vA − wA×rA)
                    real_t yphys = n.dot(vB + wBw.cross(rB) - vA - wAw.cross(rA));
                    y[cid] = yphys;
                }
                return y;
            };
            auto y_pre = contactVel(v_pre_blocks);
            auto y_np1 = contactVel(vnp1_blocks);
            auto y_pre_phys = contactVelPhys(v_pre_blocks);
            auto y_np1_phys = contactVelPhys(vnp1_blocks);
            std::cout << "[RB_DEBUG] contacts=" << contactToBlocks.size() << '\n';
            for (size_t cid = 0; cid < contactToBlocks.size() && cid < 8; ++cid) {
                auto [idxA, idxB] = contactToBlocks[cid];
                if (idxA < 0 && idxB < 0) continue; // skip static-static contacts in debug
                real_t dpre = std::abs(y_pre[cid] - y_pre_phys[cid]);
                real_t dnp1 = std::abs(y_np1[cid] - y_np1_phys[cid]);
                std::cout << "  contact " << cid
                          << ": y_pre=" << y_pre[cid]
                          << " (phys=" << y_pre_phys[cid] << ", Δ=" << dpre << ")"
                          << ", y_np1=" << y_np1[cid]
                          << " (phys=" << y_np1_phys[cid] << ", Δ=" << dnp1 << ")"
                          << ", idxA=" << idxA << ", idxB=" << idxB << '\n';
                if (idxA >= 0) {
                    int bA = blockToBody[(size_t)idxA];
                    std::cout << "    WA row (size=" << Wblocks[(size_t)idxA].cols() << ") : ";
                    for (int j = 0; j < Wblocks[(size_t)idxA].cols(); ++j) std::cout << Wblocks[(size_t)idxA](0,j) << (j+1<Wblocks[(size_t)idxA].cols()?", ":"\n");
                    VectorXr dv = vnp1_blocks[(size_t)bA] - v_pre_blocks[(size_t)bA];
                    std::cout << "    dV_A: "; for (int j=0;j<dv.size();++j) std::cout << dv[j] << (j+1<dv.size()?", ":"\n");
                }
                if (idxB >= 0) {
                    int bB = blockToBody[(size_t)idxB];
                    std::cout << "    WB row (size=" << Wblocks[(size_t)idxB].cols() << ") : ";
                    for (int j = 0; j < Wblocks[(size_t)idxB].cols(); ++j) std::cout << Wblocks[(size_t)idxB](0,j) << (j+1<Wblocks[(size_t)idxB].cols()?", ":"\n");
                    VectorXr dv = vnp1_blocks[(size_t)bB] - v_pre_blocks[(size_t)bB];
                    std::cout << "    dV_B: "; for (int j=0;j<dv.size();++j) std::cout << dv[j] << (j+1<dv.size()?", ":"\n");
                }
            }
        }
    }

    // 6) Write back final state
    m_dyn.writeStateToSystem(qnp1_blocks, vnp1_blocks);
}

}
