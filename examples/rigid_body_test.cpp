#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
#include "collision/collision_manager.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>

using namespace cardillo;

int main() {
    MPI_Init(nullptr, nullptr);
    int worldRank = 0, worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Ground plate as a static cube obstacle
    {
        PhysicsSystem::Cube ground;
        ground.center = Vector3r(0.0, 0.0, -0.5);
        ground.halfExtents = Vector3r(15.0, 15.0, 0.5); // thickness 0.5
        ground.q = Quaternion4r::Identity();
        sys.addObstacleBody(ground);
    }

    // Domino setup: a chain of upright OBB dominoes with increasing size, triggered by a moving sphere
    // Fixed parameters (no env vars)
    const int N = 36;
    const real_t h0 = (real_t)0.1;  // initial height
    const real_t dh = (real_t)1.1; // growth per domino (stronger scaling)
    const real_t aspectW = (real_t)0.45; // width ~45% of height
    const real_t aspectT = (real_t)0.10; // thickness ~10% of height (thin along chain)
    const real_t gap = (real_t)0.015;    // spacing gap between dominoes
    const real_t dens = (real_t)600.0;   // density for mass scaling

    // Build dominoes left-to-right along +X with center-to-center spacing that accounts for width growth
    real_t xCursor = -1.5; // start left of origin
    real_t prevHalfW = 0.0;
    std::vector<Vector3r> dominoCenters; dominoCenters.reserve(N);
    for (int i = 0; i < N; ++i) {
        real_t h = h0 * pow((real_t)dh, (real_t)i); // nonlinear height growth
        real_t w = aspectW * h;                   // full width (X)
        real_t t = aspectT * h;                   // full thickness (Y)
        Vector3r he(w * 0.5, t * 0.5, h * 0.5);   // half extents
        // spacing from previous: sum of half-widths + gap
        if (i == 0) {
            xCursor = -1.5 + he.x();
        } else {
            xCursor += prevHalfW + he.x() + gap;
        }
        prevHalfW = he.x();

    PhysicsSystem::Cube dom; dom.halfExtents = he;
    // Rotate 90 degrees around Z so the thin thickness (local Y) aligns with world X (chain direction)
    dom.q = Eigen::AngleAxis<real_t>((real_t)M_PI_2, Vector3r::UnitZ());
        Vector3r center(xCursor, 0.0, he.z());
        dominoCenters.push_back(center);
        // mass proportional to volume (volume = 8*he.x*he.y*he.z)
        real_t volume = (real_t)8.0 * he.x() * he.y() * he.z();
        real_t mass = std::max((real_t)0.1, dens * volume);
        sys.addRigidBody(mass, center, dom.q, Vector3r::Zero(), Vector3r::Zero(), dom);
    }

    // Rolling sphere to kick the first domino (moves along +X toward the first domino)
    {
        const real_t r = (real_t)0.06;
        const real_t vx = (real_t)2.0;
        Vector3r firstC = dominoCenters.empty() ? Vector3r(0,0,0) : dominoCenters.front();
        Vector3r sPos(firstC.x() - (real_t)0.55, 0.0, r); // start near the first domino
        Vector3r sVel(vx, 0.0, 0.0);
        const real_t m = 0.1 * dens * (4.0/3.0) * M_PI * r * r * r;
        sys.addPointMass(m, sPos, sVel, r);
    }

    // Writer (rank 0)
    std::unique_ptr<cardillo::io::VtkWriter> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriter>("vtk_out", "rigid", 50);
        writer->enableContactsOutput(true, "rigid_contacts");
    }

    // Simulate a short settling (currently only translational DOFs are integrated)
    cardillo::solver::MoreauSolver solver(sys);

    real_t t = 0.0;
    const real_t T = 5.0;
    const real_t dt = 1e-3;
    const int steps = (int)(T / dt);
    if (writer) writer->maybeWrite(0, t, sys);
    for (int k = 0; k < steps; ++k) {
        solver.stepMidpoint(dt);
        t += dt;
        if (writer) writer->maybeWrite(k+1, t, sys);
    }

    MPI_Finalize();
    return 0;
}
