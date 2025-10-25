#include "cardillo.hpp"
#include "io/vtk_writer_binary.hpp"
#include "solver/moreau.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>

using namespace cardillo;

namespace {
// Utility: compute mass from half extents and a density (with a small floor)
inline real_t massFromDensity(const Vector3r& halfExtents, real_t density, real_t minMass = (real_t)0.05) {
    const real_t volume = (real_t)8.0 * halfExtents.x() * halfExtents.y() * halfExtents.z();
    return std::max(minMass, density * volume);
}

// Ground plate as a static cube obstacle
void spawnGround(PhysicsSystem& sys, real_t halfThickness = (real_t)0.5, real_t halfSize = (real_t)15.0, real_t zTop = (real_t)0.0) {
    PhysicsSystem::Cube ground;
    ground.center = Vector3r(0.0, 0.0, zTop - halfThickness);
    ground.halfExtents = Vector3r(halfSize, halfSize, halfThickness);
    ground.q = Quaternion4r::Identity();
    sys.addObstacleBody(ground);
}

// Classic Jenga tower: layers of rectangular blocks; each layer rotated 90 degrees relative to the previous.
// Parameters:
//  - layers: number of layers in Z
//  - blockHalf: half extents of each block (in its local frame)
//  - gap: small spacing between adjacent blocks in a layer
//  - density: for mass computation
//  - baseCenter: center of the bottom layer (z will be set to blockHalf.z())
//  - blocksPerLayer: typically 3 for Jenga
//  - extraLayerGap: additional vertical gap beyond block height
void spawnJengaTower(
    PhysicsSystem& sys,
    int layers,
    const Vector3r& blockHalf,
    real_t gap,
    real_t density,
    const Vector3r& baseCenter,
    int blocksPerLayer = 3,
    real_t extraLayerGap = (real_t)0.0
) {
    if (layers <= 0 || blocksPerLayer <= 0) return;

    const real_t fullX = (real_t)2.0 * blockHalf.x(); // block length (long)
    const real_t fullY = (real_t)2.0 * blockHalf.y(); // block width (short)
    const real_t fullZ = (real_t)2.0 * blockHalf.z(); // block height

    // orientation: identity => block local x -> world x (length along X)
    // 90 deg about Z => block local x -> world y (length along Y)
    const Quaternion4r q0 = Quaternion4r::Identity();
    const Quaternion4r q90 = Quaternion4r(Eigen::AngleAxis<real_t>((real_t)M_PI_2, Vector3r::UnitZ()));

    // base z (bottom layer sits centered on baseCenter.z())
    const real_t baseZ = baseCenter.z() + blockHalf.z();

    for (int layer = 0; layer < layers; ++layer) {
        const bool alongX = (layer % 2 == 0); // even layers: long axis along X
        const Quaternion4r q = alongX ? q0 : q90;

        // We always place blocks side-by-side across the *width* (fullY).
        const real_t rowWidth = blocksPerLayer * fullY + (blocksPerLayer - 1) * gap;
        // start offset is the center of the first block relative to baseCenter
        const real_t firstOffset = -rowWidth * (real_t)0.5 + fullY * (real_t)0.5;

        const real_t z = baseZ + (real_t)layer * (fullZ + extraLayerGap);

        for (int i = 0; i < blocksPerLayer; ++i) {
            // compute block center
            Vector3r c = baseCenter;
            const real_t step = fullY + gap;
            const real_t offset = firstOffset + (real_t)i * step;

            if (alongX) {
                // block length along X, place blocks along Y (offset modifies Y)
                c = Vector3r(baseCenter.x(), baseCenter.y() + offset, z);
            } else {
                // block length along Y (rotated), place blocks along X (offset modifies X)
                c = Vector3r(baseCenter.x() + offset, baseCenter.y(), z);
            }

            PhysicsSystem::Cube blk;
            blk.halfExtents = blockHalf;
            blk.q = q;

            const real_t m = massFromDensity(blockHalf, density);
            sys.addRigidBody(m, c, blk.q, Vector3r::Zero(), Vector3r::Zero(), blk);
        }
    }
}


// Simple high-speed sphere "bullet"
void spawnBulletSphere(PhysicsSystem& sys, real_t radius, real_t density, const Vector3r& startPos, const Vector3r& velocity, real_t massScale = (real_t)1.0) {
    const real_t volume = (real_t)(4.0/3.0) * (real_t)M_PI * radius * radius * radius;
    const real_t m = massScale * std::max((real_t)0.05, density * volume);
    sys.addRigidBodySphere(m, startPos, Quaternion4r::Identity(), velocity, Vector3r::Zero(), radius);
}

// Domino tower (4-layer repeating motif), generalizable for even N x N base.
// Dominos stand on their side: long axis in-plane, thickness in-plane, height upright (Z).
void spawnDominoTowerStructure(
    PhysicsSystem& sys,
    int layers,
    int N,
    const Vector3r& half,
    real_t density,
    const Vector3r& baseCenter,
    real_t gapLong = (real_t)0.002,
    real_t extraLayerGap = (real_t)0.0
) {
    if (layers <= 0) return;
    const real_t L = (real_t)2.0 * half.x(); // long edge (in-plane)
    const real_t W = (real_t)2.0 * half.y(); // thickness (in-plane)
    const real_t H = (real_t)2.0 * half.z(); // height (upright)
    const real_t Offset = (real_t) 0.3 * W;

    // Grid spacing is domino length minus thickness, per specification
    const real_t s = std::max<real_t>((real_t)1e-6, L - (W - Offset));
    const int Ncells = std::max(2, (N / 2) * 2); // ensure even

    // A lambda function that places a domino between (i,j,k) and (i, j+1, k) or (i+1, j, k)
    auto placeDomino = [&](int i, int j, int k, bool alongY) {
        Vector3r c = baseCenter;
        // compute position offsets
        const real_t offsetX = ((real_t)i - (real_t)(Ncells - 1) * (real_t)0.5) * s - (0.5 * Offset) * alongY;
        const real_t offsetY = ((real_t)j - (real_t)(Ncells - 1) * (real_t)0.5) * s - (0.5 * Offset) * !alongY;
        c.x() += offsetX;
        c.y() += offsetY;
        const real_t z = baseCenter.z() + half.z() + (real_t)k * ( (real_t)2.0 * half.z() + extraLayerGap );
        c.z() = z;  
        PhysicsSystem::Cube domino;
        domino.halfExtents = half;
        if (alongY) {
            // rotate 90 deg about Z to align long axis along Y
            domino.q = Quaternion4r(Eigen::AngleAxis<real_t>((real_t)M_PI_2, Vector3r::UnitZ()));
            c.x() -= (L/2.0 - W/2.0);
            c.y() += (L/2.0 - W/2.0);
        } else {
            domino.q = Quaternion4r::Identity();
        }
        const real_t m = massFromDensity(half, density);
        Vector3r vel = Vector3r::Zero();
        // Give the domino in the 4th layer from the top, in the middle of y, in the most positive x a nudge
        if (i == Ncells -1 && (j == Ncells /2 || j == Ncells /2 - 1) && k == layers -4) {
            vel = Vector3r(4.0, 0.0, 1.0) * 2;
        }
        sys.addRigidBody(m, c, domino.q, vel, Vector3r::Zero(), domino);

    };

    // Place layers 
    for (int layer = 0; layer < layers; ++layer) {
        const int off = 2;
        const int k = layer;
        // Each layer is a grid of dominos
        if((layer + off) % 4 == 0) 
        {
            // place all parrallel to x from (i, j) to (i, j+1) and (i, j+2) to (i, j+3) and so on
            for (int i = 0; i < Ncells; ++i) {
                for (int j = 0; j < Ncells; j +=2) {
                    placeDomino(i, j, k, true);
                }
            }
        }else if ((layer + off) % 4 == 1)
        {
            // First place parralel along y from (i, j) to (i+1, j) and (i+2, j) to (i+3, j) and so on at bottom and top only:
            for(int i = 0; i < Ncells; i +=2)
            {
                placeDomino(i, 0, k, false);
                placeDomino(i, Ncells -1, k, false);
            }
            // Then place parralel along x along at left and right only:
            for(int j = 1; j < Ncells -1; j +=2)
            {
                placeDomino(0, j, k, true);
                placeDomino(Ncells -1, j, k, true);
            }
            // now we fill the middle parralel to y from (i+1, j+1) to (i+2, j+1) and so on
            for(int i =1; i < Ncells -1; i +=2)
            {
                for(int j =1; j < Ncells -1; j +=1)
                {
                    placeDomino(i, j, k, false);
                }
            }
        }else if ((layer + off) % 4 == 2)
        {
            // place parallel to x left and right only:
            for(int j = 0; j < Ncells; j +=2)
            {
                placeDomino(0, j, k, true);
                placeDomino(Ncells -1, j, k, true);
            }
            // place parallel to y top and bottom only:
            for(int i = 1; i < Ncells -1; i +=2)
            {
                placeDomino(i, 0, k, false);
                placeDomino(i, Ncells -1, k, false);
            }
            // now fill the middle
            for(int i =1; i < Ncells -1; i +=1)
            {
                for(int j =1; j < Ncells -1; j +=2)
                {
                    placeDomino(i, j, k, true);
                }
            }
        }else if ((layer + off) % 4 == 3)
        {
            for(int i = 0; i < Ncells; i +=2)
            {
               for(int j = 0; j < Ncells; j ++)
               {
                     placeDomino(i, j, k, false);
               }
            }
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int worldRank = 0, worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    // Load config and wire into solver and writer
    cardillo::config::Config cfg = (argc > 1)
        ? cardillo::config::ConfigReader::fromFile(argv[1])
        : cardillo::config::Config{}; // defaults from header
    
    if (argc == 0 && worldRank == 0) std::cout << "No config file provided, using defaults." << std::endl;
    
    // Initialize physics system with config (includes gravity, friction, etc.)
    PhysicsSystem sys(cfg);
    {
        const std::string exrPath = "./res/heightmaps/mountain_height.exr";
        const real_t x_dim = (real_t)50.0;  // meters across X
        const real_t y_dim = (real_t)50.0;  // meters across Y
        const real_t z_scale = (real_t)200.0; // height scaling
        const real_t min_height = (real_t)0.0; // clamp
        const Vector3r pos(0.0, 0.0, 0.0);
        const Quaternion4r q = Quaternion4r::Identity();
        sys.addObstacleHeightField(pos, q, exrPath, x_dim, y_dim, z_scale, min_height);
    }

    // // Build a domino tower near the origin (4x4 base), no bullet
    // {
    //     // Domino dims: x=length/2, y=thickness/2, z=height/2
    //     const Vector3r dominoHalf((real_t)0.048, (real_t)0.0075, (real_t)0.024); // length 8cm, thickness 1.5cm, height 4cm
    //     const real_t density = (real_t)800.0;
    //     const int layers = 7;  // 27 51
    //     const int gridN = 4;   // 16 32
    //     const Vector3r baseCenter(-3.0, 0.0, 0.0);
    //     const real_t gapLong = (real_t)0.004; // small longitudinal spacing
    //     const real_t extraLayerGap = (real_t)-0.0001;
    //     spawnDominoTowerStructure(sys, layers, gridN, dominoHalf, density, baseCenter, gapLong, extraLayerGap);
    // }

    // Drop a few rigid-body spheres onto the mountain
    {
        const real_t radius = (real_t)0.01;
        const real_t mass = (real_t)1.0;
        const int n = 5;
        for (int i = -n; i <= n; ++i) {
            for (int j = -n; j <= n; ++j) {
                Vector3r start((real_t)i * radius * 2.0, (real_t)j * radius * 2.0, (real_t)9.0);
                sys.addRigidBodySphere(mass, start, Quaternion4r::Identity(), Vector3r::Zero(), Vector3r::Zero(), radius);
            }
        }
    }

    // Setup Moreau solver
    cardillo::solver::MoreauSolver solver(sys);

    // Writer (rank 0)
    std::unique_ptr<cardillo::io::VtkWriterBinary> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriterBinary>(cfg.output_folder, cfg.output_filename_prefix, cfg.output_interval_steps);
        writer->setHeightFieldStride(cfg.output_heightfield_stride);
        if (cfg.output_write_contacts) writer->enableContactsOutput(true, cfg.output_filename_prefix + "_contacts");
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

    PetscFinalize();
    return 0;
}
