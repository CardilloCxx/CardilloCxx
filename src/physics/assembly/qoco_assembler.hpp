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

        QOCOCscMatrix toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A);
        QOCOCscMatrix toQocoCSC(const Eigen::VectorX<real_t>& v);
        QOCOFloat* toQocoVector(const Eigen::VectorX<real_t>& v);

        QOCOCscMatrix P();
        QOCOCscMatrix A();
        QOCOCscMatrix G();

        QOCOFloat* c();
        QOCOFloat* b();
        QOCOFloat* h();

    private:
        cardillo::physics::DynamicsAssembler* m_dyn{nullptr};
};

} // namespace cardillo::physics::assembly
