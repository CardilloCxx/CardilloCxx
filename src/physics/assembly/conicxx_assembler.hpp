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

// Standard form of ConicXX (as Clarabel): min 1/2 x'Px + q'x + 1/2 z_E'C z_E s.t. Ax - C z_E + s = b,
// s in K = Zero x Nonneg x SOC x ... x SOC, with C the compliance of the zero-cone rows.
//
// Compliant rows (default): x = u, P = M, the spring and damper rows W_g', W_gamma' are zero-cone
// rows with compliances C = (Ĉ_g, Ĉ_gamma); the bilateral percussions are pi = -z_E.
// Augmented (compliant_rows = false): x = (u, pi_g, pi_gamma), P = blkdiag(M, Ĉ_g, Ĉ_gamma),
// rows [W_g' Ĉ_g 0; W_gamma' 0 Ĉ_gamma] with C = 0. Both have the same solution u, z.
class ConicxxAssembler {
   public:
    ConicxxAssembler(physics::DynamicsAssembler& dyn, bool compliant_rows) : m_dyn(&dyn), m_compliant_rows(compliant_rows) {}

    bool compliantRows() const { return m_compliant_rows; }

    const SparseMatrix<Eigen::ColMajor>& P(real_t dt, real_t theta);
    VectorXr& q(real_t dt, real_t theta);

    const SparseMatrix<Eigen::ColMajor>& A(real_t dt, real_t theta);
    VectorXr& b(real_t dt, real_t theta);

    // Compliance of the zero-cone rows (all zero in the augmented form).
    const VectorXr& compliance(real_t dt, real_t theta);

    const conicxx::ConeSpec& coneSpec();

    VectorXr computeSmu() const;

   private:
    // Diagonal (Ĉ_g, Ĉ_gamma) of the spring and damper rows.
    VectorXr bilateralCompliance(real_t dt, real_t theta) const;

    physics::DynamicsAssembler* m_dyn{nullptr};
    bool m_compliant_rows{true};

    SparseMatrix<Eigen::ColMajor> m_P_cache;
    SparseMatrix<Eigen::ColMajor> m_A_cache;
    Eigen::VectorX<real_t> m_q_cache;
    Eigen::VectorX<real_t> m_b_cache;
    Eigen::VectorX<real_t> m_compliance_cache;
    conicxx::ConeSpec m_cone_spec_cache;
};

}  // namespace cardillo::physics::assembly
