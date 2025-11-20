#include "cardillo.hpp"
#include "io/vtk_writer_binary.hpp"
#include "solver/moreau.hpp"
#include "solver/dual_stoermer_verlet.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>
#include <chrono>
#include "misc/progress/ProgressBar.hpp"
#include <csignal>
#include <cstdlib>
#include <iostream>

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
#include "scenes/wilberforce/WilberforcePendulum.hpp"
#include "scenes/hangbride/HangbrideScene.hpp"
#include "scenes/softbody/SoftbodyTestScene.hpp"
#include "scenes/parcel/ParcelScene.hpp"
#include "scenes/rodAssembly/RodAssemblyScene.hpp"
#include "scenes/discreteRod/DiscreteRodScene.hpp"
#include "scenes/gears/GearsScene.hpp"
#include "scenes/pendulum/PendulumScene.hpp"

using namespace cardillo;

static PhysicsSystem sys(cardillo::config::Config{}); 

void printTimingsAtExit(int sig) {
    sys.timings().printBreakdown(std::cout);
    std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int worldRank = 0, worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    std::signal(SIGINT, printTimingsAtExit);
    Eigen::setNbThreads(1);

    // Load config and wire into solver and writer
    cardillo::config::Config cfg = (argc > 1)
        ? cardillo::config::ConfigReader::fromFile(argv[1])
        : cardillo::config::Config{}; // defaults from header
    
    if (argc == 0 && worldRank == 0) std::cout << "No config file provided, using defaults." << std::endl;
    
    sys.setConfig(cfg);

    // HeightmapScene scene;
    // DominoScene scene;
    // SpringTestScene scene;
    // RodAssemblyScene scene;
    // NetScene scene;
    // WilberforcePendulumScene scene;
    // HangbrideScene scene;
    // SoftbodyTestScene scene;
    // ChainScene scene;
    // PainleveScene scene;
    // RotatingBallScene scene;
    // RailScene scene;
    // DzhanibekovScene scene;
    // ParcelScene scene;
    // DiscreteRodScene scene;
    // ConstraintTestScene scene;
    // GearsScene scene;
    PendulumScene scene;

    scene.populate(sys);

    // Setup Moreau solver
    cardillo::solver::MoreauSolver solver(sys);
    // cardillo::solver::DualStoermerVerletSolver solver(sys);

    // Writer (rank 0)
    std::unique_ptr<cardillo::io::VtkWriterBinary> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriterBinary>(cfg.output_folder, cfg.output_filename_prefix, cfg.output_interval_steps);
        writer->setHeightFieldStride(cfg.output_heightfield_stride);
        if (cfg.output_write_contacts) writer->enableContactsOutput(true, cfg.output_filename_prefix + "_contacts");
        writer->enableSpringsOutput(true, cfg.output_filename_prefix + std::string("_springs"));
    }

    real_t t = 0.0;
    const real_t T = cfg.sim_T;
    const real_t dt = cfg.sim_dt;
    const int steps = (int)(T / dt);
    if (writer) writer->maybeWrite(0, t, sys);

    auto t0 = std::chrono::steady_clock::now();
    std::unique_ptr<cardillo::misc::ProgressBar> pbar;

    if (worldRank == 0) {
        pbar = std::make_unique<cardillo::misc::ProgressBar>(static_cast<std::size_t>(steps), std::cout);
        pbar->set_description("Simulating");
    }
    {
        auto totalScope = sys.timings().scope(cardillo::misc::TimingManager::TimerId::Total);
        for (int k = 0; k < steps; ++k) {
            scene.updateScene(sys, t, dt);
            solver.stepMidpoint(dt);
            t += dt;
            if (writer) writer->maybeWrite(k+1, t, sys);
            if (pbar) pbar->update(1);
        }
    }
    if (pbar) pbar->close();
    if (worldRank == 0) {
        sys.timings().printBreakdown(std::cout);
    }

    PetscFinalize();
    return 0;
}