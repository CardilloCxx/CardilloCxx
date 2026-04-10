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
    PjAssembler(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}
    ~PjAssembler();

    bool buildAndFactorS(real_t dt, real_t theta, bool implicitGyro = false, bool lambdaTheta = false);
    VectorXr rhs(real_t dt, real_t theta) const;
    VectorXr solveS(const VectorXr& rhs_ext) const;
    const TripletMatrix& S() const { return m_S; }

   private:
    TripletMatrix m_S;
    std::optional<Eigen::SparseLU<CscMatrix>> m_S_sparse_lu;
    cardillo::physics::DynamicsAssembler* m_dyn{nullptr};
    cardillo::config::Config m_cfg;
};

}  // namespace cardillo::physics::assembly

#endif