#include "cardillo.hpp"
#include "physics/dynamics_assembler.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace cardillo;

int main() {
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Two masses: 1.0 and 2.0
    sys.addPointMass(1.0, Vector3r(0,0,0), Vector3r(0,0,0), 0.05);
    sys.addPointMass(2.0, Vector3r(0,0,0), Vector3r(0,0,0), 0.05);

    // Mass inverse blocks should be diag(1,1,1) and diag(1/2,1/2,1/2)
    cardillo::physics::DynamicsAssembler dyn(sys);
    dyn.assignDofs();
    dyn.refreshState();
    const auto& offV = dyn.bodyVelOffsets();
    const auto& MinvDiag = dyn.MinvDiag();
    assert((int)offV.size() == 3); // 2 bodies + 1
    // Body 0
    for (int i = 0; i < 3; ++i) assert(std::abs(MinvDiag[offV[0] + i] - 1.0) < 1e-12);
    // Body 1
    for (int i = 0; i < 3; ++i) assert(std::abs(MinvDiag[offV[1] + i] - 0.5) < 1e-12);

    // Force blocks should be [m*g; m*g]
    const auto& f = dyn.fVec();
    Vector3r g(0,0,-9.81);
    for (int i = 0; i < 3; ++i) assert(std::abs(f[offV[0] + i] - 1.0 * g[i]) < 1e-12);
    for (int i = 0; i < 3; ++i) assert(std::abs(f[offV[1] + i] - 2.0 * g[i]) < 1e-12);

    std::cout << "matrix_tests OK\n";
    return 0;
}
