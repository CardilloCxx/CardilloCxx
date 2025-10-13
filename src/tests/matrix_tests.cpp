#include "cardillo.hpp"
#include "physics/dynamics_assembler.hpp"
#include <iostream>
#include <cassert>

using namespace cardillo;

// Tiny utility to get dense matrix from sparse
static MatrixXXr dense(const Eigen::SparseMatrix<real_t>& S) {
    return MatrixXXr(S);
}

int main() {
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Two masses: 1.0 and 2.0
    sys.addPointMass(1.0, Vector3r(0,0,0), Vector3r(0,0,0), 0.05);
    sys.addPointMass(2.0, Vector3r(0,0,0), Vector3r(0,0,0), 0.05);

    // Mass matrix should be diag(1,1,1, 2,2,2)
    cardillo::physics::DynamicsAssembler dyn(sys);
    auto M = dyn.M();
    auto Md = dense(M);
    assert(Md.rows() == 6 && Md.cols() == 6);
    for (int i = 0; i < 3; ++i) assert(std::abs(Md(i,i) - 1.0) < 1e-12);
    for (int i = 3; i < 6; ++i) assert(std::abs(Md(i,i) - 2.0) < 1e-12);

    // Derive M diagonal from dyn.M()
    VectorXr Mdiag(Md.rows());
    for (int i = 0; i < Md.rows(); ++i) Mdiag[i] = Md(i,i);
    assert(Mdiag.size() == 6);
    for (int i = 0; i < 6; ++i) assert(std::abs(Mdiag[i] - Md(i,i)) < 1e-12);

    // Force vector should be [m*g; m*g]
    VectorXr f = dyn.f();
    assert(f.size() == 6);
    for (int i = 0; i < 3; ++i) assert(std::abs(f[i] - 1.0 * Vector3r(0,0,-9.81)[i]) < 1e-12);
    for (int i = 0; i < 3; ++i) assert(std::abs(f[3+i] - 2.0 * Vector3r(0,0,-9.81)[i]) < 1e-12);

    std::cout << "matrix_tests OK\n";
    return 0;
}
