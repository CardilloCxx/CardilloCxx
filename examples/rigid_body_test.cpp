#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
#include "collision/collision_manager.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>

using namespace cardillo;

int main(int argc, char** argv) {
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

    // // Block laying on the ground
    // {
    //     PhysicsSystem::Cube block;
    //     block.center = Vector3r(0.0, 0.0, 0.24);
    //     block.halfExtents = Vector3r(1.0, 1.0, 0.25); // thickness 0.25
    //     Eigen::AngleAxis<real_t> rotAngleAxis( (real_t)(0.5 * M_PI / 180.0), Vector3r::UnitX() );
    //     block.q = Quaternion4r(rotAngleAxis);
    //     sys.addObstacleBody(block);
    // }

    // 2D pyramid stack of flat rectangular boxes (lying on the ground), with spacing
    // Layout (top to bottom), centered around x=0:
    //             [####]
    //        [####]   [####]
    //    [####]   [####]   [####]
    // [####]   [####]   [####]   [####]
//     {
//         // Dimensions (full extents) of each box
//         const real_t boxLenX = (real_t)0.6;  // length along X
//         const real_t boxLenY = (real_t)0.2;  // depth along Y (into page)
//         const real_t boxLenZ = (real_t)0.1;  // thickness along Z (height)
//         const Vector3r half(boxLenX*(real_t)0.5, boxLenY*(real_t)0.5, boxLenZ*(real_t)0.5);
// 
//         // Spacing between boxes in the same row (gap along X)
//         const real_t gapX = (real_t)0.3;
//         // Vertical spacing between rows (additional to box height)
//         const real_t gapZ = (real_t)0.01;
//         // Base row count
//         const int baseCount = 10;
//         // Density for mass, or choose a fixed mass
//         const real_t density = (real_t)60.0;
// 
//         auto boxMass = [&](const Vector3r& he){
//             const real_t volume = (real_t)8.0 * he.x() * he.y() * he.z(); // full extents = 2*he
//             return std::max((real_t)0.05, density * volume);
//         };
// 
//         // Compute base row starting Z so that bottom boxes rest on ground (ground top at z=0)
//         const real_t baseZ = half.z();
// 
//         // Build rows from bottom (row=0 has baseCount) to top (row=baseCount-1 has 1)
//         for (int row = 0; row < baseCount; ++row) {
//             int count = baseCount - row;
// 
//             // Total width of the row including gaps between boxes
//             const real_t rowWidth = (real_t)count * boxLenX + (real_t)(count - 1) * gapX;
//             // Leftmost X so that the row is centered around x=0
//             const real_t xLeft = -(real_t)0.5 * rowWidth + (real_t)0.5 * boxLenX;
// 
//             // Z position for this row
//             const real_t z = baseZ + (real_t)row * (boxLenZ + gapZ);
// 
//             for (int i = 0; i < count; ++i) {
//                 const real_t x = xLeft + (real_t)i * (boxLenX + gapX);
//                 const real_t y = (real_t)0.0; // centered along Y
// 
//                 PhysicsSystem::Cube cube;
//                 cube.halfExtents = half; // flat rectangle
//                 cube.q = Quaternion4r::Identity(); // lying flat on ground
// 
//                 const Vector3r center(x, y, z);
//                 const real_t m = boxMass(half);
//                 sys.addRigidBody(m, center, cube.q, Vector3r::Zero(), Vector3r::Zero(), cube);
//             }
//         }
//     }

//     // Domino setup: a chain of upright OBB dominoes with increasing size, triggered by a moving sphere
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



    // Load config and wire into solver and writer
    cardillo::config::Config cfg = (argc > 1)
        ? cardillo::config::ConfigReader::fromFile(argv[1])
        : cardillo::config::Config{}; // defaults from header
    cardillo::solver::MoreauSolver solver(sys, cfg);

    // Writer (rank 0)
    std::unique_ptr<cardillo::io::VtkWriter> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriter>("vtk_out", "rigid", cfg.output_interval_steps);
        writer->enableContactsOutput(true, "rigid_contacts");
    }

    real_t t = 0.0;
    const real_t T = cfg.sim_T;
    const real_t dt = cfg.sim_dt;
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
