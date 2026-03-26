#pragma once

#include "misc/types.hpp"
#include "physics/world.hpp"
#include "physics/assembly/dynamics_assembler.hpp"
#include "solver/solver_base.hpp"
#include "solver/projected_jacobi.hpp"
#include <memory>

namespace cardillo::integration {

class IntegrationBase {
public:
    explicit IntegrationBase(cardillo::World& world, cardillo::physics::DynamicsAssembler& dyn,
                             cardillo::misc::TimingManager& timings,
                             cardillo::solver::WarmstartProvider* warmstart)
        : m_world(world), m_dyn(dyn), m_timings(timings) {
        // Create default numeric solver based on config (currently ProjectedJacobi)
        m_solver = std::make_unique<cardillo::solver::ProjectedJacobiSolver>(m_dyn, world.config(),
                                                                             warmstart);
    }

    virtual ~IntegrationBase() = default;
    virtual void step(real_t dt) = 0;

    // Access to the underlying numeric solver (may be null)
    cardillo::solver::SolverBase* solver() { return m_solver.get(); }
    const cardillo::solver::SolverBase* solver() const { return m_solver.get(); }

    void explicitPositionUpdate(cardillo::World& world, real_t h);
    void linearImplicitPositionUpdate(cardillo::World& world, real_t h);

protected:
    cardillo::World& m_world;
    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::misc::TimingManager& m_timings;
    std::unique_ptr<cardillo::solver::SolverBase> m_solver;
};

} // namespace cardillo::integration
