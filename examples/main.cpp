#include "cardillo.hpp"
#include "io/vtk_writer_binary.hpp"
#include "solver/solver_base.hpp"
#include "solver/moreau.hpp"
#include "solver/dual_stoermer_verlet.hpp"
#include <iostream>
#include <iomanip>
#include <Eigen/Geometry>
#include <chrono>
#include "misc/progress/ProgressBar.hpp"
#include <csignal>
#include <cstdlib>
#include <iostream>

// Scenes
#include "scenes/SceneBase.hpp"
#include "scenes/heightmap/HeightmapScene.hpp"
#include "scenes/domino/DominoScene.hpp"
#include "scenes/cardhouse/CardhouseScene.hpp"
#include "scenes/painleve/painleveScene.hpp"
#include "scenes/rotating_ball/RotatingBallScene.hpp"
#include "scenes/chain/ChainScene.hpp"
#include "scenes/rail/RailScene.hpp"
#include "scenes/dzhanibekov/dzhanibekov.hpp"
#include "scenes/springTest/SpringTestScene.hpp"
#include "scenes/rodAssembly/RodAssemblyScene.hpp"
#include "scenes/net/NetScene.hpp"
#include "scenes/euler_disk/EulerDiskScene.hpp"
#include "scenes/wilberforce/WilberforcePendulum.hpp"
#include "scenes/hangbride/HangbrideScene.hpp"
#include "scenes/softbody/SoftbodyTestScene.hpp"
#include "scenes/parcel/ParcelScene.hpp"
#include "scenes/rodAssembly/RodAssemblyScene.hpp"
#include "scenes/discreteRod/DiscreteRodScene.hpp"
#include "scenes/slinky/SlinkyScene.hpp"
#include "scenes/unilateral/UnilateralScene.hpp"
#include "scenes/woodpecker/WoodpeckerScene.hpp"
#include "scenes/spaghetti/SpaghettiScene.hpp"
#include "scenes/metronome/MetronomeScene.hpp"
#include "scenes/stacked_spheres/StackedSpheresScene.hpp"
#include "scenes/strandbeest/StrandbeestScene.hpp"
// #include "scenes/gears/GearsScene.hpp"         
// #include "scenes/pendulum/PendulumScene.hpp"  
#include "scenes/leaningTower/LeaningTowerScene.hpp"
#include "scenes/constraintTest/ConstraintTestScene.hpp"
#include "scenes/cantilever/CantileverScene.hpp"
#include "scenes/fabric/FabricScene.hpp"
#include "scenes/heightmap/HeightmapScene.hpp"
#include "scenes/double_pendulum/DoublePendulumScene.hpp"

using namespace cardillo;

static World sys(cardillo::config::Config{}); 

void printTimingsAtExit(int sig) {
    sys.timings().printBreakdown(std::cout);
    std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    const int worldRank = 0;

    std::signal(SIGINT, printTimingsAtExit);
    Eigen::setNbThreads(1);

    // Load config and wire into solver and writer
    cardillo::config::Config cfg = (argc > 1)
        ? cardillo::config::ConfigReader::fromFile(argv[1])
        : cardillo::config::Config{}; // defaults from header
    
    if (argc <= 1 && worldRank == 0) std::cout << "No config file provided, using defaults." << std::endl;
    
    sys.setConfig(cfg);

    // Construct all available scenes and select the one matching cfg.scene_name
    std::vector<std::unique_ptr<SceneBase>> scenes;
    scenes.emplace_back(std::make_unique<HeightmapScene>());
    scenes.emplace_back(std::make_unique<StackedSpheresScene>());
    scenes.emplace_back(std::make_unique<DominoScene>());
    scenes.emplace_back(std::make_unique<CardhouseScene>());
    scenes.emplace_back(std::make_unique<SpringTestScene>());
    scenes.emplace_back(std::make_unique<RodAssemblyScene>());
    scenes.emplace_back(std::make_unique<NetScene>());
    scenes.emplace_back(std::make_unique<WilberforcePendulumScene>());
    scenes.emplace_back(std::make_unique<HangbrideScene>());
    scenes.emplace_back(std::make_unique<SoftbodyTestScene>());
    scenes.emplace_back(std::make_unique<ChainScene>());
    scenes.emplace_back(std::make_unique<PainleveScene>());
    scenes.emplace_back(std::make_unique<RotatingBallScene>());
    scenes.emplace_back(std::make_unique<RailScene>());
    scenes.emplace_back(std::make_unique<DzhanibekovScene>());
    scenes.emplace_back(std::make_unique<ParcelScene>());
    scenes.emplace_back(std::make_unique<DiscreteRodScene>());
    scenes.emplace_back(std::make_unique<SlinkyScene>());
    scenes.emplace_back(std::make_unique<EulerDiskScene>());
    scenes.emplace_back(std::make_unique<UnilateralScene>());
    scenes.emplace_back(std::make_unique<WoodpeckerScene>());
    scenes.emplace_back(std::make_unique<SpaghettiScene>());
    scenes.emplace_back(std::make_unique<MetronomeScene>());
    scenes.emplace_back(std::make_unique<StrandbeestScene>());
    // scenes.emplace_back(std::make_unique<GearsScene>());   
    scenes.emplace_back(std::make_unique<LeaningTowerScene>());
    // scenes.emplace_back(std::make_unique<PendulumScene>());  
    scenes.emplace_back(std::make_unique<ConstraintTestScene>());
    scenes.emplace_back(std::make_unique<CantileverScene>());
    scenes.emplace_back(std::make_unique<FabricScene>());
    scenes.emplace_back(std::make_unique<HeightmapScene>());
    scenes.emplace_back(std::make_unique<DoublePendulumScene>());

    SceneBase* selected = nullptr;
    for (auto& s : scenes) {
        if (cfg.scene_name == s->sceneName()) {
            selected = s.get();
            break;
        }
    }
    if (!selected) {
        if (worldRank == 0) 
        {
            if (cfg.scene_name == "none-specified") {
                std::cerr << "No scene_name specified in config. Available scenes are:" << std::endl;
                for (auto& s : scenes) {
                    std::cerr << "  - " << s->sceneName() << std::endl;
                }
            }
            std::cerr << "Unknown scene_name '" << cfg.scene_name << "'" << std::endl;
        }
        PetscFinalize();
        return EXIT_FAILURE;
    }
    std::cout << "Selected scene: " << selected->sceneName() << std::endl;
    SceneBase& scene = *selected;
    cardillo::physics::PhysicsEngine engine(sys);
    cfg.output_filename_prefix = scene.sceneName();
    scene.populate(engine);

    // Setup solver based on config
    std::unique_ptr<cardillo::solver::SolverBase> solver;
    if (cfg.solver == cardillo::config::SolverType::StoermerVerlet) {
        if (worldRank == 0) {
            std::cout << "[Warning] Dual Stoermer-Verlet is deprecated; prefer Moreau (solver.name=moreau).\n";
        }
        solver = std::make_unique<cardillo::solver::DualStoermerVerletSolver>(sys);
    } else {
        solver = std::make_unique<cardillo::solver::MoreauSolver>(sys, cfg.moreau_theta);
    }

    std::cout << "[Info] Selected solver: " 
              << ((cfg.solver == cardillo::config::SolverType::StoermerVerlet) ? "Dual Stoermer-Verlet" : "Moreau") 
              << std::endl;

    if (cfg.solver == cardillo::config::SolverType::Moreau) {
        std::cout << "[Info] Moreau settings: theta=" << cfg.moreau_theta
                  << ", implicit_gyroscopy=" << (cfg.moreau_implicit_gyroscopy ? "true" : "false")
                  << std::endl;
    }
    
    // Writer (rank 0)
    std::unique_ptr<cardillo::io::VtkWriterBinary> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriterBinary>(cfg.output_folder, cfg.output_filename_prefix, cfg.output_interval_steps);
        writer->setHeightFieldStride(cfg.output_heightfield_stride);
        if (cfg.output_write_contacts) writer->enableContactsOutput(true, cfg.output_filename_prefix + std::string("_contacts"));
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
            scene.updateScene(engine, t, dt);
            solver->stepMidpoint(dt);
            if (worldRank == 0) {
                sys.writeTrackedStateToCsv(t + dt);
            }
            t += dt;
            if (writer) writer->maybeWrite(k+1, t, sys);
            if (pbar) {
                int jorIters = solver->lastProjectedJacobiIterations();
                if (jorIters >= 0) pbar->set_postfix("jor=" + std::to_string(jorIters) + "        ");
                pbar->update(1);
            }
        }
    }
    if (pbar) pbar->close();
    if (worldRank == 0) {
        sys.timings().printBreakdown(std::cout);
    }

    PetscFinalize();
    return 0;
}