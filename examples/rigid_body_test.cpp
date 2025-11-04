#include "cardillo.hpp"
#include "io/vtk_writer_binary.hpp"
#include "solver/moreau.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>

// Scenes
#include "scenes/heightmap/HeightmapScene.hpp"
#include "scenes/domino/DominoScene.hpp"
#include "scenes/painleve/painleveScene.hpp"
#include "scenes/rotating_ball/RotatingBallScene.hpp"
#include "scenes/chain/ChainScene.hpp"
#include "scenes/rail/RailScene.hpp"
#include "scenes/dzhanibekov/dzhanibekov.hpp"
#include "scenes/dzhanibekov/dzhanibekov.hpp"
#include "scenes/springTest/SpringTestScene.hpp"
#include "scenes/rodAssembly/RodAssemblyScene.hpp"
#include "scenes/net/NetScene.hpp"
#include "scenes/hangbride/HangbrideScene.hpp"
#include "scenes/softbody/SoftbodyTestScene.hpp"
#include "scenes/parcel/ParcelScene.hpp"

using namespace cardillo;

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
    // HeightmapScene scene;
    // DominoScene scene;
    // SpringTestScene scene;
    // RodAssemblyScene scene;
    // NetScene scene;
    // HangbrideScene scene;
    // SoftbodyTestScene scene;
    // ChainScene scene;
    // PainleveScene scene;
    // RotatingBallScene scene;
    // RailScene scene;
    // DzhanibekovScene scene;
    ParcelScene scene;

    scene.populate(sys);

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
        scene.maybeSpawnParcel(sys, t);
        scene.resetPosition(sys);
        // sys.markForcesDirty();
        // sys.markStructureDirty();
        // sys.markStateDirty();
        // sys.
        if (writer) writer->maybeWrite(k+1, t, sys);
    }

    if (worldRank == 0) {
        std::cout << "Simulation finished — exiting cleanly." << std::endl;
    }
    PetscFinalize();
    return 0;
}
