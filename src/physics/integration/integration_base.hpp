#pragma once

#include <memory>
#include "../solver/projected_jacobi.hpp"
#include "../solver/solver_base.hpp"
#include "misc/types.hpp"
#include "physics/assembly/dynamics_assembler.hpp"
#include "physics/world.hpp"

namespace cardillo::integration {

class IntegrationBase {
   public:
    explicit IntegrationBase(cardillo::World& world, cardillo::solver::SolverBase& solver, cardillo::physics::DynamicsAssembler& dyn, cardillo::misc::TimingManager& timings,
                             cardillo::config::Config& config)
        : m_world(world), m_solver(solver), m_dyn(dyn), m_timings(timings), m_config(config) {}

    virtual ~IntegrationBase() = default;
    virtual void step(real_t dt) = 0;

    void explicitPositionUpdate(cardillo::World& world, real_t h);
    void linearImplicitPositionUpdate(cardillo::World& world, real_t h);

   protected:
    cardillo::World& m_world;
    cardillo::physics::DynamicsAssembler& m_dyn;
    cardillo::misc::TimingManager& m_timings;
    cardillo::solver::SolverBase& m_solver;
    cardillo::config::Config& m_config;
};

}  // namespace cardillo::integration
