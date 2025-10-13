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
    real_t pradius = 0.05;
    sys.addPointMass(mass, pos, vel, pradius);
    }

    // Mass matrix is now provided by DynamicsAssembler via MoreauSolver if needed

    // Configure visual planes (e.g., z=0 ground)
    PhysicsSystem::Plane ground;
    ground.center = Vector3r(0,0,-1);
    ground.normal = Vector3r(0,0,1);
    ground.up = Vector3r(0,1,0);
    ground.sizeX = 5.0; ground.sizeY = 5.0;
    sys.addRigidBody(ground);

    // Add a cube visual off to the side
    PhysicsSystem::Cube cube;
    cube.center = Vector3r(2.0, 0.0, 0.0);
    cube.halfExtents = Vector3r(0.5, 0.5, 0.5);
    sys.addRigidBody(cube);

    // Step a few times and export VTK every 5 steps
    cardillo::io::VtkWriter writer("vtk_out", "points", 5);
    real_t dt = 0.01;
    cardillo::solver::MoreauSolver solver(sys);
    for (int k = 0; k < 100; ++k) {
        solver.stepMidpoint(dt);
        writer.maybeWrite(k, sys);

        // Print state: grab the first visual point via ECS
        if (k % 20 == 0) {
            const auto& reg = sys.ecs();
            auto view = reg.view<PhysicsSystem::C_VisualObject,
                                 PhysicsSystem::C_PointVisualTag,
                                 PhysicsSystem::C_Position3,
                                 PhysicsSystem::C_LinearVelocity3>();
            for (auto [e, pos, vel] : view.each()) {
                (void)e;
                std::cout << "Step " << k << ": sample x[0] = " << pos.value.transpose()
                          << ", v[0] = " << vel.value.transpose() << "\n";
                break;
            }
        }
    }

    return 0;
}
