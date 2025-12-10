#pragma once

#include "../misc/types.hpp"

namespace cardillo::solver {

class SolverBase {
public:
    virtual ~SolverBase() = default;
    virtual void stepMidpoint(real_t dt) = 0;
};

} // namespace cardillo::solver
