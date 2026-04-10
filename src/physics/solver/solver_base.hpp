#pragma once

#include <Eigen/Dense>
#include "../../misc/types.hpp"

namespace cardillo {
namespace solver {

class SolverBase {
   public:
    virtual ~SolverBase() = default;

    // Solve the system in-place or return a result vector. Signature mirrors ProjectedJacobi.
    virtual VectorXr solve(real_t dt, real_t theta) = 0;

    // Diagnostics
    virtual int lastIterations() const { return m_last_iters; }
    virtual real_t lastError() const { return (real_t)0; }

    virtual const char* name() const { return "SolverBase"; }

   protected:
    int m_last_iters{-1};
};

}  // namespace solver
}  // namespace cardillo
