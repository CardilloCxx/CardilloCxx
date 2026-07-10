#pragma once

#include "io/vtk_writer.hpp"
#include "physics/integration/integration_base.hpp"
#include "physics/integration/moreau.hpp"

#include <Eigen/Geometry>

#include <csignal>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "config/config.hpp"
#include "physics/ecs_types.hpp"
#include "scenes/SceneBase.hpp"

namespace cardillo::examples {

namespace detail {
inline cardillo::physics::PhysicsEngine* g_engine = nullptr;

inline void printTimingsAtExit(int sig) {
    (void)sig;
    if (g_engine) g_engine->timings().printBreakdown(std::cout);
    std::exit(EXIT_FAILURE);
}
}  // namespace detail

template <typename SceneType>
int runExample(int argc, char** argv) {
    std::signal(SIGINT, detail::printTimingsAtExit);
    Eigen::setNbThreads(1);

    cardillo::config::Config cfg = (argc > 1) ? cardillo::config::ConfigReader::fromFile(argv[1]) : cardillo::config::Config{};

    if (argc <= 1) std::cout << "No config file provided, using defaults." << std::endl;

    SceneType scene;
    cfg.scene_name = scene.sceneName();
    cfg.output_filename_prefix = scene.sceneName();

    cardillo::physics::PhysicsEngine engine(cfg);
    detail::g_engine = &engine;
    scene.populate(engine);

    real_t t = 0.0;
    const real_t dt = cfg.sim_dt;
    while (!engine.isFinished()) {
        scene.updateScene(engine, t, dt);
        engine.step();
        t += dt;
    }
    engine.timings().printBreakdown(std::cout);

    if (std::getenv("CARDILLO_DUMP_STATE")) {
        real_t totalKE = 0;
        real_t posNormSum = 0;
        real_t velNormSum = 0;
        int numBodies = 0;
        const auto& reg = engine.ecs();
        auto view = reg.view<cardillo::C_BodyIndex>();
        for (auto e : view) {
            ++numBodies;
            totalKE += engine.getKineticEnergy(e);
            posNormSum += engine.getPosition(e).norm();
        }
        std::cerr << std::setprecision(15) << "[STATE-DUMP] t=" << t << " numBodies=" << numBodies << " totalKE=" << totalKE << " posNormSum=" << posNormSum << std::endl;
    }

    return 0;
}

}  // namespace cardillo::examples

#define CARDILLO_DEFINE_EXAMPLE_MAIN(SceneType) \
    int main(int argc, char** argv) { return ::cardillo::examples::runExample<SceneType>(argc, argv); }
