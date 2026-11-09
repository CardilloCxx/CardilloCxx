#include "moreau.hpp"
#include <iostream>
#include "../solver/contact_residual.hpp"
#include "../solver/projected_jacobi.hpp"

namespace cardillo::integration {

void MoreauIntegrator::step(real_t dt) {
    m_dyn.refreshState();

    // First drift to the intermediate configuration q_{n+theta} at time t_n + (1 - theta) dt.
    // Prescribed constraint rates and kinematic drivers (in updateStateDependentTerms()) are
    // advanced to the same time, so that g and W are evaluated consistently at t_{n+theta}.
    explicitPositionUpdate(m_world, (1.0 - m_config.moreau_theta) * dt);
    m_dyn.advancePrescribedMotion((1.0 - m_config.moreau_theta) * dt);
    m_dyn.updateStateDependentTerms(dt);
    const VectorXr u = m_solver.solve(dt, m_config.moreau_theta);
    if (m_config.pj_final_check == config::PjFinalCheck::NaturalMap) {
        const real_t error = solver::contactNaturalMapResidual(m_dyn, u, m_config.pj_tol_abs, m_config.pj_tol_rel);
        if (m_config.debug_pj) std::cout << "Final natural-map error: " << error << std::endl;
    }
    m_dyn.writeVelocityToSystem(u, dt);
    explicitPositionUpdate(m_world, m_config.moreau_theta * dt);
    m_dyn.advancePrescribedMotion(m_config.moreau_theta * dt);
}

}  // namespace cardillo::integration