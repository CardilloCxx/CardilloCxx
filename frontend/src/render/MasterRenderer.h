#pragma once

#include <filament/Engine.h>
#include <filament/Renderer.h>
#include <filament/View.h>
#include <filament/Camera.h>
#include <filament/SwapChain.h>
#include <filament/ColorGrading.h>
#include <utils/Entity.h>
#include <memory>
#include "world/GameWorld.h"
#include "world/EntityPicker.h"
#include "camera/CameraController.h"

namespace frontend::render {

class MasterRenderer {
public:
    MasterRenderer(const platform::InputState* inputState);
    ~MasterRenderer();

    bool initialize(void* nativeWindowHandle, int width, int height);
    void resize(int width, int height);
    
    // Pass the world in so the renderer can swap scenes if needed
    bool render(GameWorld* world, double elapsedSeconds);

    filament::Engine* getEngine() const { return engine_; }
    EntityPicker* getEntityPicker() const { return entityPicker_.get(); }

private:
    const platform::InputState* inputState_ = nullptr;
    
    filament::Engine* engine_ = nullptr;
    filament::Renderer* renderer_ = nullptr;
    filament::SwapChain* swapChain_ = nullptr;
    
    // Main View
    filament::View* mainView_ = nullptr;
    filament::Camera* camera_ = nullptr;
    utils::Entity cameraEntity_;
    filament::ColorGrading* colorGrading_ = nullptr;
    filament::Renderer::ClearOptions mainClearOptions_{};

    // Subsystems
    std::unique_ptr<camera::CameraController> cameraController_;
    std::unique_ptr<EntityPicker> entityPicker_;
    
    int width_ = 1280;
    int height_ = 720;
    void* nativeWindowHandle_ = nullptr;

    uint32_t lastHoveredEntityId = 0;

};

} // namespace filament_example::render