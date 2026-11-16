#include "moreau.hpp"
#include <chrono>
#include "../solver/natural_map.hpp"
#include "../synchronization/derived_entity_sync.hpp"

namespace cardillo::integration {

void MoreauIntegrator::step(real_t dt) {
    m_dyn.refreshState();

    // First drift to the intermediate configuration q_{n+theta} at time t_n + (1 - theta) dt.
    // Prescribed constraint rates and kinematic drivers (in updateStateDependentTerms()) are
    // advanced to the same time, so that g and W are evaluated consistently at t_{n+theta}.
    explicitPositionUpdate(m_world, (1.0 - m_config.moreau_theta) * dt);
    m_dyn.advancePrescribedMotion((1.0 - m_config.moreau_theta) * dt);
    m_dyn.updateStateDependentTerms(dt);
    const auto solveStart = std::chrono::steady_clock::now();
    const VectorXr u = m_solver.solve(dt, m_config.moreau_theta);
    m_lastSolveSeconds = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - solveStart).count();
    if (m_config.solver_final_check == config::FinalCheck::NaturalMap) m_lastResidual = solver::stepResidual(m_dyn, u, m_solver.contactLaw(), dt, m_config.moreau_theta);
    m_dyn.writeVelocityToSystem(u, dt);
    explicitPositionUpdate(m_world, m_config.moreau_theta * dt);
    m_dyn.advancePrescribedMotion(m_config.moreau_theta * dt);
    // The cached RigidState (read by outputs, applyForceAt, new constraints) now holds q_{n+1}, u_{n+1}.
    physics::DerivedEntitySync::refreshRigidStates(m_world);
}

}  // namespace cardillo::integration