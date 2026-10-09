#include "moreau.hpp"
#include "../solver/projected_jacobi.hpp"

namespace cardillo::integration {

void MoreauIntegrator::step(real_t dt) {
    m_dyn.refreshState();
    // Position-level force law evaluated at the beginning of the step: g(t_n, q_n) before the first drift.

    // First drift to the intermediate configuration q_{n+theta} at time t_n + (1 - theta) dt.
    // Prescribed constraint rates and kinematic drivers (in updateStateDependentTerms()) are
    // advanced to the same time, so that g and W are evaluated consistently at t_{n+theta}.
    explicitPositionUpdate(m_world, (1.0 - m_config.moreau_theta) * dt);
    m_dyn.advancePrescribedMotion((1.0 - m_config.moreau_theta) * dt);
    m_dyn.updateStateDependentTerms(dt);
    m_dyn.writeVelocityToSystem(m_solver.solve(dt, m_config.moreau_theta), dt);
    explicitPositionUpdate(m_world, m_config.moreau_theta * dt);
    m_dyn.advancePrescribedMotion(m_config.moreau_theta * dt);
}

}  // namespace cardillo::integration