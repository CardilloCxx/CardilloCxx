#ifndef PJ_ASSEMBLER_H
#define PJ_ASSEMBLER_H

#pragma once

#include <optional>
#include "../../config/config.hpp"
#include "../../misc/triplet_matrix.hpp"
#include "../../misc/types.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {
class PjAssembler

{
   public:
    PjAssembler(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}
    ~PjAssembler();

    bool buildAndFactorS(real_t dt, real_t theta, bool implicitGyro = false, bool lambdaTheta = false);
    VectorXr rhs(real_t dt, real_t theta) const;
    VectorXr solveS(const VectorXr& rhs_ext) const;
    const TripletMatrix& S() const { return m_S; }

   private:
    TripletMatrix m_S;
    std::optional<Eigen::SparseLU<CscMatrix>> m_S_sparse_lu;
    physics::DynamicsAssembler* m_dyn{nullptr};
    config::Config m_cfg;

    // Structural sparsity-pattern cache: SparseLU::analyzePattern() (fill-reducing column ordering +
    // symbolic elimination) is far more expensive than factorize() and only needs to be redone when
    // S's nonzero *pattern* changes, not merely its values (which change every step via dt/state).
    // Mirrors CondensedAssembler's bilateral-order cache (see condensed_assembler.cpp).
    bool m_hasCachedPattern{false};
    bool m_cachedIsSymmetric{false};
    Eigen::Index m_cachedRows{0}, m_cachedCols{0};
    std::vector<CscMatrix::StorageIndex> m_cachedOuterIndex;
    std::vector<CscMatrix::StorageIndex> m_cachedInnerIndex;
};

}  // namespace cardillo::physics::assembly

#endif