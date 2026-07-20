#pragma once

#include <vector>

#include <Eigen/SparseCore>

#include <conicxx/cone_spec.h>
#include "../../misc/triplet_matrix.hpp"
#include "../../misc/types.hpp"

namespace cardillo {
class World;
}  // namespace cardillo
namespace cardillo::physics {
class DynamicsAssembler;
}  // namespace cardillo::physics

namespace cardillo::physics::assembly {

// Builds the same standard-form problem as ClarabelAssembler
// (min 1/2 x'Px + q'x s.t. Ax + s = b, s in K = Zero x Nonneg x SOC x ... x SOC)
// since ConicXX and Clarabel share that exact convention.
class ConicxxAssembler {
   public:
    explicit ConicxxAssembler(physics::DynamicsAssembler& dyn) : m_dyn(&dyn) {}

    const SparseMatrix<Eigen::ColMajor>& P(real_t dt, real_t theta);
    VectorXr& q(real_t dt, real_t theta);

    const SparseMatrix<Eigen::ColMajor>& A(real_t dt, real_t theta);
    VectorXr& b(real_t dt, real_t theta);

    const conicxx::ConeSpec& coneSpec();

    VectorXr computeSmu() const;

   private:
    physics::DynamicsAssembler* m_dyn{nullptr};

    SparseMatrix<Eigen::ColMajor> m_P_cache;
    SparseMatrix<Eigen::ColMajor> m_A_cache;
    Eigen::VectorX<real_t> m_q_cache;
    Eigen::VectorX<real_t> m_b_cache;
    conicxx::ConeSpec m_cone_spec_cache;
};

}  // namespace cardillo::physics::assembly
