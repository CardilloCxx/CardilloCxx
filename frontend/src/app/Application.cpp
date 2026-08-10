#include "app/Application.h"
#include "world/GameWorld.h"
#include "render/MasterRenderer.h"
#include "physics/CardilloPhysicsBridge.h"
#include "scene/Scene.h"

#include <chrono>
#include <iostream>
#include <thread>

namespace frontend::app {

Application::Application()
    : platform_(std::make_unique<platform::Platform>(platform::WindowConfig{})) {}

int Application::run() {
    if (!platform_->initialize()) {
        std::cerr << "Unable to initialize the platform window\n";
        return 1;
    }

    auto renderer = std::make_unique<render::MasterRenderer>(&platform_->inputState());
    if (!renderer->initialize(platform_->nativeWindowHandle(), platform_->width(), platform_->height())) {
        std::cerr << "Unable to initialize the Filament renderer\n";
        return 2;
    }

    auto world = std::make_unique<GameWorld>(*renderer->getEngine(), *renderer->getEntityPicker());
    world->initialize();

    physicsBridge_ = std::make_unique<frontend::physics::CardilloPhysicsBridge>(*world);
    physicsBridge_->initialize("examples/scenes/jenga/scene.config");
    // frontend::scene::setupCornellBox(*world);

    const double targetFPS = 60.0;
    const std::chrono::duration<double> frameTargetDuration(1.0 / targetFPS);
    constexpr double maxComputationTimePerFrame = 0.015;

    auto lastFrame = std::chrono::steady_clock::now();
 
    while (!platform_->shouldClose()) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = now - lastFrame;

        if (elapsed < frameTargetDuration) {
            std::this_thread::yield(); 
            continue; 
        }

        const double elapsedSeconds = std::chrono::duration<double>(elapsed).count();
        lastFrame = now;

        platform_->pollEvents();
        physicsBridge_->step(elapsedSeconds, maxComputationTimePerFrame);
        world->update(elapsedSeconds);
        renderer->render(world.get(), elapsedSeconds);
    }

    world.reset(); 
    renderer.reset();

    return 0;
}

} // namespace frontend::app