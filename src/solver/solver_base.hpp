#pragma once

#include "../misc/types.hpp"
#include <Eigen/Dense>

namespace cardillo { namespace solver {

class SolverBase {
public:
    virtual ~SolverBase() = default;

    // Solve the system in-place or return a result vector. Signature mirrors ProjectedJacobi.
    virtual VectorXr solve(VectorXr& rhs, real_t tol = (real_t)1e-5) = 0;

    // Diagnostics
    virtual int lastIterations() const { return -1; }
    virtual real_t lastError() const { return (real_t)0; }

    virtual const char* name() const { return "SolverBase"; }
};

}} // namespace cardillo::solver
