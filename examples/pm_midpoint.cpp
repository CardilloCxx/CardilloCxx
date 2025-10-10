#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
#include <iostream>

using namespace cardillo;

int main() {
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Create a ring of point masses with tangential initial velocities
    const int N = 80;
    real_t radius = 1.0;
    real_t mass = 1.0;
    real_t v0 = 1.0; // tangential speed
    for (int i = 0; i < N; ++i) {
        real_t theta = 2.0 * M_PI * i / N;
        Vector3r pos(radius * std::cos(theta), radius * std::sin(theta), 0.0);
        // tangential direction in XY plane
        Vector3r tan(-std::sin(theta), std::cos(theta), 0.0);
        Vector3r vel = v0 * tan;
        sys.addPointMass(mass, pos, vel);
    }

    auto M = sys.assembleMassMatrix();
    std::cout << "Mass matrix M (dense) size = " << M.rows() << " x " << M.cols() << "\n\n";

    // Step a few times and export VTK each step
    cardillo::io::VtkWriter writer("vtk_out", "points", 1);
    real_t dt = 0.01;
    for (int k = 0; k < 100; ++k) {
        cardillo::solver::midpointStep(sys, dt);
        writer.maybeWrite(k, sys);

        // Print state
        if (k % 20 == 0) {
            const auto& masses = sys.masses();
            std::cout << "Step " << k << ": sample x[0] = " << masses[0].x.transpose()
                      << ", v[0] = " << masses[0].v.transpose() << "\n";
        }
    }

    return 0;
}
