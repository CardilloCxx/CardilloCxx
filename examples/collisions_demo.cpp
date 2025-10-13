#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
// collision/collision.hpp not needed directly when using VtkWriter::maybeWriteAll

using namespace cardillo;

int main() {
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Ground plane (z=0)
    PhysicsSystem::Plane ground;
    ground.center = Vector3r(0,0,0);
    ground.normal = Vector3r(0,0,1);
    ground.up = Vector3r(0,1,0);
    ground.sizeX = 4.0; ground.sizeY = 4.0;
    sys.addRigidBody(ground);

    // Tilted finite plane near y=-1
    PhysicsSystem::Plane tilt;
    tilt.center = Vector3r(0,-1.0,0.0);
    tilt.normal = Vector3r(0, std::sqrt(0.5), std::sqrt(0.5));
    tilt.up = Vector3r(1,0,0);
    tilt.sizeX = 1.0; tilt.sizeY = 1.0;
    sys.addRigidBody(tilt);

    // Rotated cube obstacle
    PhysicsSystem::Cube cube;
    cube.center = Vector3r(0.5, 0.0, 0.2);
    cube.halfExtents = Vector3r(0.3, 0.2, 0.2);
    real_t ang = M_PI / 6.0;
    Matrix33r Rz; Rz = Eigen::AngleAxis<real_t>(ang, Vector3r::UnitZ());
    cube.R = Rz;
    sys.addRigidBody(cube);

    // Two spheres heading towards each other
    real_t m = 1.0;
    real_t r = 0.25;
    sys.addPointMass(m, Vector3r(1.0, 0.0, 0.5), Vector3r(-1.0, 0.0, 0.0), r);
    sys.addPointMass(m, Vector3r(-1.0, 0.0, 0.5), Vector3r(1.0, 0.0, 0.0), r);

    
    sys.addPointMass(m, Vector3r(0.0, 0.0, 1.0), Vector3r(0.0, 0.0, 1.0), r);

    // Writers
    cardillo::io::VtkWriter writer("vtk_out", "scene", 1);
    writer.enableContactsOutput(true, "contacts");

    // Simulate
    const int steps = 200;
    const real_t dt = 0.01;

    cardillo::solver::MoreauSolver solver(sys);
    for (int k = 0; k < steps; ++k) {
        solver.stepMidpoint(dt);
        writer.maybeWrite(k, sys);
    }

    return 0;
}
