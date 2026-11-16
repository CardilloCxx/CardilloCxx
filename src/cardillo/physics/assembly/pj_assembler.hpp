#ifndef PJ_ASSEMBLER_H
#define PJ_ASSEMBLER_H

#pragma once

#include <optional>
#include "../../config/config.hpp"
#include "../../misc/fixed_pattern_matrix.hpp"
#include "../../misc/types.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {
class PjAssembler

{
   public:
    PjAssembler(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}
    ~PjAssembler();

    bool buildAndFactorS(real_t dt, real_t theta, bool implicitGyro = false);
    VectorXr rhs(real_t dt, real_t theta) const;
    VectorXr solveS(const VectorXr& rhs_ext) const;
    // Solves S X = B for all columns of B at once (one factorization, many right-hand sides).
    MatrixXXr solveSColumns(const MatrixXXr& B) const;

   private:
    // Saddle-point matrix S with a fixed pattern (see buildAndFactorS()) and its factorization;
    // m_analyzedSymmetric is the symmetry flag of the last analyzePattern().
    misc::FixedPatternMatrix m_S;
    std::optional<Eigen::SparseLU<CscMatrix>> m_S_sparse_lu;
    bool m_analyzedSymmetric{false};
    physics::DynamicsAssembler* m_dyn{nullptr};
    config::Config m_cfg;
};

}  // namespace cardillo::physics::assembly

#endif