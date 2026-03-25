#pragma once

#include "misc/types.hpp"
#include "physics/world.hpp"
#include "physics/assembly/dynamics_assembler.hpp"
#include "solver/projected_jacobi.hpp"

namespace cardillo::integration {

class IntegrationBase {
public:
    explicit IntegrationBase(cardillo::World& world)
        : m_world(world), m_dyn(world),
          m_pj(m_dyn, world.config(), world.config().pj_warmstart ? world.warmstartProvider() : nullptr) {}

    virtual ~IntegrationBase() = default;
    virtual void step(real_t dt) = 0;

    // Returns last Projected-Jacobi iteration count if available, else -1
    virtual int lastProjectedJacobiIterations() const { return -1; }

    static void explicitPositionUpdate(cardillo::World& world, real_t h);
    static void linearImplicitPositionUpdate(cardillo::World& world, real_t h);

protected:
    cardillo::World& m_world;
    cardillo::physics::DynamicsAssembler m_dyn;
    cardillo::solver::ProjectedJacobiSolver m_pj;
};

} // namespace cardillo::integration
