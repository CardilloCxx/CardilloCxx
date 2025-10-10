#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../misc/types.hpp"
#include "../physics/physics_system.hpp"

namespace cardillo::solver {

// Midpoint rule for unconstrained translation-only point masses
void midpointStep(cardillo::PhysicsSystem& sys, real_t dt);

}

