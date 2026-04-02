#include "io/vtk_writer_binary.hpp"
#include "physics/integration/integration_base.hpp"
#include "physics/integration/moreau.hpp"
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

static cardillo::physics::PhysicsEngine* g_engine = nullptr;

void printTimingsAtExit(int sig) {
    if (g_engine) g_engine->timings().printBreakdown(std::cout);
    std::exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    std::signal(SIGINT, printTimingsAtExit);
    Eigen::setNbThreads(1);

    // Load config and wire into solver and writer
    cardillo::config::Config cfg = (argc > 1)
        ? cardillo::config::ConfigReader::fromFile(argv[1])
        : cardillo::config::Config{}; // defaults from header

    // print working drectory
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
    
    if (argc <= 1) std::cout << "No config file provided, using defaults." << std::endl;

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
        if (cfg.scene_name == "none-specified") {
            std::cerr << "No scene_name specified in config. Available scenes are:" << std::endl;
            for (auto& s : scenes) {
                std::cerr << "  - " << s->sceneName() << std::endl;
            }
        }
        std::cerr << "Unknown scene_name '" << cfg.scene_name << "'" << std::endl;

        PetscFinalize();
        return EXIT_FAILURE;
    }
    
    std::cout << "Selected scene: " << selected->sceneName() << std::endl;

    // Construct engine from config
    SceneBase& scene = *selected;
    cfg.output_filename_prefix = scene.sceneName();
    cardillo::physics::PhysicsEngine engine(cfg);
    g_engine = &engine;
    scene.populate(engine);

    real_t t = 0.0;
    const real_t dt = cfg.sim_dt;
    while (!engine.isFinished()) {
        scene.updateScene(engine, t, dt);
        engine.step();
        t += dt;
    }
    engine.timings().printBreakdown(std::cout);

    PetscFinalize();
    return 0;
}