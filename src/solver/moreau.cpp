#include "moreau.hpp"
#include "projected_jacobi.hpp"
#include <vector>
#include <algorithm>
#include <limits>

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

    // 2) Midpoint position and preliminary velocity (external forces only)
    std::vector<VectorXr> q_mid_blocks(Nb), v_pre_blocks(Nb);
    for (int b = 0; b < Nb; ++b) {
        q_mid_blocks[(size_t)b] = qnBlocks[(size_t)b] + (real_t)0.5 * dt * vnBlocks[(size_t)b];
        v_pre_blocks[(size_t)b] = vnBlocks[(size_t)b] + dt * a_ext_blocks[(size_t)b];
    }

    // 3) Evaluate contacts at midpoint positions: write q_mid into ECS (v stays vn)
    m_dyn.writeStateToSystem(q_mid_blocks, vnBlocks);
    m_dyn.refreshCollisions();

    // 4) Solve normal impulses p using W(q_mid), G, and preliminary velocity v_pre
    std::vector<VectorXr> vnp1_blocks = m_pj.iterateWithPreliminaryVelocity(v_pre_blocks, (real_t)1e-5);

    // 5) Final position: q_{n+1} = q_mid + (h/2) * u_{n+1}
    std::vector<VectorXr> qnp1_blocks(Nb);
    for (int b = 0; b < Nb; ++b) qnp1_blocks[(size_t)b] = q_mid_blocks[(size_t)b] + (real_t)0.5 * dt * vnp1_blocks[(size_t)b];

    // 6) Write back final state
    m_dyn.writeStateToSystem(qnp1_blocks, vnp1_blocks);
}

}
