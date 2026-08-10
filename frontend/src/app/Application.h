#pragma once

#include <memory>

#include "platform/Platform.h"
#include "render/MasterRenderer.h"
#include "scene/Scene.h"
#include "physics/CardilloPhysicsBridge.h"

namespace frontend::app {

class Application {
public:
    Application();
    int run();

private:
    void processInput(double elapsedSeconds);

    std::unique_ptr<platform::Platform> platform_;
    std::unique_ptr<render::MasterRenderer> renderer_;
    std::unique_ptr<frontend::physics::CardilloPhysicsBridge> physicsBridge_;
};

} // namespace frontend::app
