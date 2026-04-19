#pragma once

#include <vector>

#include <Eigen/SparseCore>

#include <clarabel.hpp>
#include "../../misc/triplet_matrix.hpp"
#include "../../misc/types.hpp"

namespace cardillo {
class World;
}  // namespace cardillo
namespace cardillo::physics {
class DynamicsAssembler;
}  // namespace cardillo::physics

namespace cardillo::physics::assembly {

class ClarabelAssembler {
   public:
    explicit ClarabelAssembler(cardillo::physics::DynamicsAssembler& dyn) : m_dyn(&dyn) {}

    const SparseMatrix<Eigen::ColMajor>& P(real_t dt, real_t theta);
    VectorXr& q(real_t dt, real_t theta);

    const SparseMatrix<Eigen::ColMajor>& A(real_t dt, real_t theta);
    VectorXr& b(real_t dt, real_t theta);

    const std::vector<clarabel::SupportedConeT<double>>& cones();
    std::size_t coneCount() const { return m_cones_cache.size(); }

    VectorXr computeSmu() const;

   private:
    cardillo::physics::DynamicsAssembler* m_dyn{nullptr};

    SparseMatrix<Eigen::ColMajor> m_P_cache;
    SparseMatrix<Eigen::ColMajor> m_A_cache;
    Eigen::VectorX<real_t> m_q_cache;
    Eigen::VectorX<real_t> m_b_cache;
    std::vector<clarabel::SupportedConeT<double>> m_cones_cache;
};

}  // namespace cardillo::physics::assembly
