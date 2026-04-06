#include "moreau.hpp"
#include "../solver/projected_jacobi.hpp"

namespace cardillo::integration {

void MoreauIntegrator::step(real_t dt)
{
    m_dyn.refreshState();

    explicitPositionUpdate(m_world, (1.0 - m_config.moreau_theta) * dt);
    m_dyn.writeVelocityToSystem(m_solver.solve(dt, m_config.moreau_theta), dt);
    explicitPositionUpdate(m_world, m_config.moreau_theta * dt);
}

}