#pragma once

#include <vector>
#include <Eigen/SparseCore>
#include "../../misc/types.hpp"
#include "../../misc/triplet_matrix.hpp"
#include <qoco/qoco_linalg.h>
#include <qoco/structs.h>

namespace cardillo::physics {
class World;
class DynamicsAssembler;
}

namespace cardillo::physics::assembly {

class QocoAssembler {

    public:
        QocoAssembler(cardillo::physics::DynamicsAssembler& dyn) : m_dyn(&dyn) {}

        QOCOCscMatrix* toQocoCSC(SparseMatrix<Eigen::ColMajor>& A, QOCOCscMatrix& view);
        QOCOFloat* toQocoVector(Eigen::VectorX<real_t>& v);

        QOCOCscMatrix* P(real_t dt, real_t theta);
        QOCOCscMatrix* A(real_t dt, real_t theta);
        QOCOCscMatrix* G(real_t dt, real_t theta);

        QOCOFloat* c(real_t dt, real_t theta);
        QOCOFloat* b(real_t dt, real_t theta);
        QOCOFloat* h(real_t dt, real_t theta);

        VectorXr computeSmu();

    private:
        cardillo::physics::DynamicsAssembler* m_dyn{nullptr};

        // Keep all sparse/vector storage alive while QOCO consumes raw pointers.
        SparseMatrix<Eigen::ColMajor> m_A_cache;
        SparseMatrix<Eigen::ColMajor> m_G_cache;
        SparseMatrix<Eigen::ColMajor> m_P_cache;

        Eigen::VectorX<real_t> m_c_cache;
        Eigen::VectorX<real_t> m_b_cache;
        Eigen::VectorX<real_t> m_h_cache;

        QOCOCscMatrix m_P_view{};
        QOCOCscMatrix m_A_view{};
        QOCOCscMatrix m_G_view{};
};

} // namespace cardillo::physics::assembly
