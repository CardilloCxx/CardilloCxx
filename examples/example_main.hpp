#pragma once

#include "io/vtk_writer_binary.hpp"
#include "physics/integration/integration_base.hpp"
#include "physics/integration/moreau.hpp"

#include <Eigen/Geometry>

#include <csignal>
#include <cstdlib>
#include <iostream>

#include "config/config.hpp"
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

    return 0;
}

}  // namespace cardillo::examples

#define CARDILLO_DEFINE_EXAMPLE_MAIN(SceneType) \
    int main(int argc, char** argv) { return ::cardillo::examples::runExample<SceneType>(argc, argv); }
