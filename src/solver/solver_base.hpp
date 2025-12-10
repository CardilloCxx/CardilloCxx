#pragma once

#include "../misc/types.hpp"

namespace cardillo::solver {

class SolverBase {
public:
    virtual ~SolverBase() = default;
    virtual void stepMidpoint(real_t dt) = 0;
    // Returns last Projected-Jacobi iteration count if available, else -1
    virtual int lastProjectedJacobiIterations() const { return -1; }
};

} // namespace cardillo::solver
