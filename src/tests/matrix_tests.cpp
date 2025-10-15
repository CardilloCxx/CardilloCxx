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
    dyn.loadStateFromSystem();
    dyn.refreshState();
    const auto& MinvBlocks = dyn.MinvBlocks();
    assert((int)MinvBlocks.size() == 2);
    for (int i = 0; i < 3; ++i) assert(std::abs(MinvBlocks[0](i,i) - 1.0) < 1e-12);
    for (int i = 0; i < 3; ++i) assert(std::abs(MinvBlocks[1](i,i) - 0.5) < 1e-12);

    // Force blocks should be [m*g; m*g]
    const auto& fBlocks = dyn.f();
    assert((int)fBlocks.size() == 2);
    Vector3r g(0,0,-9.81);
    for (int i = 0; i < 3; ++i) assert(std::abs(fBlocks[0][i] - 1.0 * g[i]) < 1e-12);
    for (int i = 0; i < 3; ++i) assert(std::abs(fBlocks[1][i] - 2.0 * g[i]) < 1e-12);

    std::cout << "matrix_tests OK\n";
    return 0;
}
