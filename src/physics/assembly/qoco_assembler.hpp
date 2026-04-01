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
        QOCOCscMatrix toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A);
        QOCOCscMatrix toQocoCSC(const Eigen::VectorX<real_t>& v);
        QOCOVectorf* toQocoVectorf(const Eigen::VectorX<real_t>& v);

        QOCOCscMatrix P();
        QOCOCscMatrix A();
        QOCOCscMatrix G();

        QOCOVectorf* c();
        QOCOVectorf* b();
        QOCOVectorf* h();

    private:
        cardillo::physics::World* m_world{nullptr};
        cardillo::physics::DynamicsAssembler* m_dyn{nullptr};
};

} // namespace cardillo::physics::assembly
