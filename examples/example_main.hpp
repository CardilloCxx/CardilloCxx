#pragma once

#include "io/vtk_writer.hpp"
#include "physics/integration/integration_base.hpp"
#include "physics/integration/moreau.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "config/config.hpp"
#include "physics/ecs_types.hpp"
#include "scenes/SceneBase.hpp"

namespace cardillo::examples {

namespace detail {
inline physics::PhysicsEngine* g_engine = nullptr;

// Largest violation |g| over all perfect constraint rows (rows with zero compliance), evaluated at
// the current state. Diagnostic for CARDILLO_DUMP_STATE; compliant rows are excluded since their
// elongation is physical.
inline real_t maxPerfectConstraintViolation(physics::PhysicsEngine& engine) {
    // Constraints read the cached RigidState, which still holds the intermediate configuration of
    // the last step; refresh it so that the violation is evaluated at the current configuration.
    auto& reg = engine.ecs();
    for (auto e : reg.view<const C_Position3>()) RigidBody::updateState(reg, e);
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

    config::Config cfg = (argc > 1) ? config::ConfigReader::fromFile(argv[1]) : config::Config{};

    if (argc <= 1) std::cout << "No config file provided, using defaults." << std::endl;

    SceneType scene;
    cfg.scene_name = scene.sceneName();
    cfg.output_filename_prefix = scene.sceneName();

    physics::PhysicsEngine engine(cfg);
    detail::g_engine = &engine;
    scene.populate(engine);

    real_t t = 0.0;
    const real_t dt = cfg.sim_dt;
    const bool dumpState = std::getenv("CARDILLO_DUMP_STATE") != nullptr;
    real_t maxViolationOverRun = 0;
    while (!engine.isFinished()) {
        // Call through the deprecated overload so that scenes still overriding updateScene(engine,
        // t, dt) keep working; its default forwards to updateScene(engine, t).
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        static_cast<SceneBase&>(scene).updateScene(engine, t, dt);
#pragma GCC diagnostic pop
        engine.step();
        t += dt;
        if (dumpState) maxViolationOverRun = std::max(maxViolationOverRun, detail::maxPerfectConstraintViolation(engine));
    }
    engine.timings().printBreakdown(std::cout);

    if (std::getenv("CARDILLO_DUMP_STATE")) {
        real_t totalKE = 0;
        real_t posNormSum = 0;
        real_t velNormSum = 0;
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

    return 0;
}

}  // namespace cardillo::examples

#define CARDILLO_DEFINE_EXAMPLE_MAIN(SceneType) \
    int main(int argc, char** argv) { return ::cardillo::examples::runExample<SceneType>(argc, argv); }
