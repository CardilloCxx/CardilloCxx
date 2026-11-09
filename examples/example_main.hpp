#pragma once

#include <cardillo/physics/integration/integration_base.hpp>
#include <cardillo/physics/integration/moreau.hpp>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <cardillo/config/config.hpp>
#include <cardillo/physics/ecs_types.hpp>
#include "scenes/SceneBase.hpp"

namespace cardillo::examples {

namespace detail {
// Set by SIGINT; the main loop then stops and prints the timings (the handler itself may only
// touch this flag).
inline volatile std::sig_atomic_t g_interrupted = 0;

// Largest violation |g| over all perfect constraint rows (rows with zero compliance), evaluated at
// the current state. Diagnostic for CARDILLO_DUMP_STATE; compliant rows are excluded since their
// elongation is physical.
inline real_t maxPerfectConstraintViolation(physics::PhysicsEngine& engine) {
    real_t maxViolation = 0;
    for (const auto& pattern : engine.world().constraintPatterns()) {
        if (!pattern) continue;
        const auto res = pattern->getConstraint();
        for (int i = 0; i < (int)res.Crows.size() && i < (int)res.positionError.size(); ++i) {
            if (res.Crows[i] == (real_t)0) maxViolation = std::max(maxViolation, std::abs(res.positionError[i]));
        }
    }
    return maxViolation;
}

inline void onInterrupt(int /*sig*/) { g_interrupted = 1; }
}  // namespace detail

template <typename SceneType>
int runExample(int argc, char** argv) {
    std::signal(SIGINT, detail::onInterrupt);
    Eigen::setNbThreads(1);

    config::Config cfg = (argc > 1) ? config::ConfigReader::fromFile(argv[1]) : config::Config{};

    if (argc <= 1) std::cout << "No config file provided, using defaults." << std::endl;

    SceneType scene;
    // Config values win, so that config variants of one scene can write to different files.
    if (cfg.scene_name == config::Config{}.scene_name) cfg.scene_name = scene.sceneName();
    if (cfg.output_filename_prefix.empty()) cfg.output_filename_prefix = scene.sceneName();

    physics::PhysicsEngine engine(cfg);
    scene.populate(engine);

    real_t t = 0.0;
    const real_t dt = cfg.sim_dt;
    const bool dumpState = std::getenv("CARDILLO_DUMP_STATE") != nullptr;
    real_t maxViolationOverRun = 0;
    while (!engine.isFinished() && !detail::g_interrupted) {
        scene.updateScene(engine, t);
        engine.step();
        t += dt;
        if (dumpState) maxViolationOverRun = std::max(maxViolationOverRun, detail::maxPerfectConstraintViolation(engine));
    }
    if (detail::g_interrupted) std::cout << "\nInterrupted at t = " << t << std::endl;
    engine.timings().printBreakdown(std::cout);

    if (std::getenv("CARDILLO_DUMP_STATE")) {
        real_t totalKE = 0;
        real_t posNormSum = 0;
        int numBodies = 0;
        const auto& reg = engine.ecs();
        auto view = reg.view<C_BodyIndex>();
        for (auto e : view) {
            ++numBodies;
            totalKE += engine.getKineticEnergy(e);
            posNormSum += engine.getPosition(e).norm();
        }
        std::cerr << std::setprecision(15) << "[STATE-DUMP] t=" << t << " numBodies=" << numBodies << " totalKE=" << totalKE << " posNormSum=" << posNormSum << " maxConstraintViolation=" << maxViolationOverRun
                  << " finalConstraintViolation=" << detail::maxPerfectConstraintViolation(engine) << std::endl;
    }

    return detail::g_interrupted ? EXIT_FAILURE : 0;
}

}  // namespace cardillo::examples

#define CARDILLO_DEFINE_EXAMPLE_MAIN(SceneType) \
    int main(int argc, char** argv) { return ::cardillo::examples::runExample<SceneType>(argc, argv); }
