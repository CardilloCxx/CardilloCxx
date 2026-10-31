#pragma once

#include <vector>
#include "../../config/config.hpp"
#include "../../misc/block_diagonal.hpp"
#include "../../misc/types.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

class PgsAssembler {
   public:
    PgsAssembler(physics::DynamicsAssembler& dyn, const config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}

    // Convenience overload: computes ufree() internally. Prefer the explicit overload below when
    // the caller already needs ufree() itself, to avoid computing it twice per solve.
    VectorXr rhs(real_t dt, real_t theta) const { return rhs(dt, theta, ufree(dt, theta)); }
    VectorXr rhs(real_t dt, real_t theta, const VectorXr& u_free) const;
    VectorXr ufree(real_t dt, real_t theta) const;
    BlockDiagonal Dinv(real_t dt, real_t theta) const;
    BlockDiagonal DinvDiag(real_t dt, real_t theta) const;

    Eigen::SparseMatrix<real_t> W() const;
    VectorXr C(real_t dt, real_t theta) const;

   private:
    physics::DynamicsAssembler* m_dyn{nullptr};
    config::Config m_cfg;
};

}  // namespace cardillo::physics::assembly
