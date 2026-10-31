#pragma once

#include <vector>

#include <Eigen/SparseCore>

#include "../../../misc/triplet_matrix.hpp"
#include "../../../misc/types.hpp"

namespace cardillo {
class World;
}  // namespace cardillo
namespace cardillo::physics {
class DynamicsAssembler;
}  // namespace cardillo::physics

namespace cardillo::physics::assembly {

// Cone program min ½xᵀPx + qᵀx s.t. Ax + s = b, s ∈ K shared by the Clarabel and SCS backends.
// Row order of A/b: zero cone (springs, dampers), nonnegative cone (frictionless contacts), then
// one 3-dimensional second-order cone per frictional contact.
class ConeProgramAssembler {
   public:
    explicit ConeProgramAssembler(physics::DynamicsAssembler& dyn) : m_dyn(&dyn) {}

    const SparseMatrix<Eigen::ColMajor>& P(real_t dt, real_t theta);
    VectorXr& q(real_t dt, real_t theta);

    const SparseMatrix<Eigen::ColMajor>& A(real_t dt, real_t theta);
    VectorXr& b(real_t dt, real_t theta);

    VectorXr computeSmu() const;

   private:
    physics::DynamicsAssembler* m_dyn{nullptr};

    SparseMatrix<Eigen::ColMajor> m_P_cache;
    SparseMatrix<Eigen::ColMajor> m_A_cache;
    Eigen::VectorX<real_t> m_q_cache;
    Eigen::VectorX<real_t> m_b_cache;
};

}  // namespace cardillo::physics::assembly
