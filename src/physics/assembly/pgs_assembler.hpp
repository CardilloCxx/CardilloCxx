#pragma once

#include <vector>
#include "../../config/config.hpp"
#include "../../misc/block_diagonal.hpp"
#include "../../misc/types.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

class PgsAssembler {
   public:
    PgsAssembler(cardillo::physics::DynamicsAssembler& dyn, const cardillo::config::Config& cfg) : m_dyn(&dyn), m_cfg(cfg) {}

    VectorXr rhs(real_t dt, real_t theta) const;
    VectorXr ufree(real_t dt, real_t theta) const;
    BlockDiagonal Dinv(real_t dt, real_t theta) const;
    BlockDiagonal DinvDiag(real_t dt, real_t theta) const;

    Eigen::SparseMatrix<real_t> W() const;
    VectorXr C(real_t dt, real_t theta) const;

   private:
    cardillo::physics::DynamicsAssembler* m_dyn{nullptr};
    cardillo::config::Config m_cfg;
};

}  // namespace cardillo::physics::assembly
