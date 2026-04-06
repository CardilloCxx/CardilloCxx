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

        QOCOCscMatrix* toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A);
        QOCOCscMatrix* toQocoCSC(const Eigen::VectorX<real_t>& v);
        QOCOFloat* toQocoVector(const Eigen::VectorX<real_t>& v);

        QOCOCscMatrix* P(real_t dt, real_t theta);
        QOCOCscMatrix* A(real_t dt, real_t theta);
        QOCOCscMatrix* G(real_t dt, real_t theta);

        QOCOFloat* c(real_t dt, real_t theta);
        QOCOFloat* b(real_t dt, real_t theta);
        QOCOFloat* h(real_t dt, real_t theta);

        VectorXr computeSmu();

    private:
        cardillo::physics::DynamicsAssembler* m_dyn{nullptr};

        // Keep sparse storage alive while QOCO reads CSC pointers (zero-copy path).
        SparseMatrix<Eigen::ColMajor> m_A_cache;
        SparseMatrix<Eigen::ColMajor> m_G_cache;
};

} // namespace cardillo::physics::assembly
