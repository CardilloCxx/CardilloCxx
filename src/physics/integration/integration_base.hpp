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
    explicit IntegrationBase(World& world, solver::SolverBase& solver, physics::DynamicsAssembler& dyn, misc::TimingManager& timings,
                             config::Config& config)
        : m_world(world), m_dyn(dyn), m_timings(timings), m_solver(solver), m_config(config) {}

    virtual ~IntegrationBase() = default;
    virtual void step(real_t dt) = 0;

    void explicitPositionUpdate(World& world, real_t h);
    void linearImplicitPositionUpdate(World& world, real_t h);

   protected:
    World& m_world;
    physics::DynamicsAssembler& m_dyn;
    misc::TimingManager& m_timings;
    solver::SolverBase& m_solver;
    config::Config& m_config;
};

}  // namespace cardillo::integration
