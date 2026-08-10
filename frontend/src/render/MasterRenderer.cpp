#include "MasterRenderer.h"
#include <filament/Viewport.h>
#include <filament/ColorGrading.h>
#include <filament/Renderer.h>
#include <filament/Engine.h>
#include <filament/Camera.h>
#include <filament/View.h>
#include <utils/EntityManager.h>
#include <iostream>
#include <cmath>

namespace frontend::render {

MasterRenderer::MasterRenderer(const platform::InputState* inputState) 
    : inputState_(inputState) {}

MasterRenderer::~MasterRenderer() {
    if (!engine_) return;

    if (colorGrading_) engine_->destroy(colorGrading_);
    entityPicker_.reset(); // Destroy custom systems before engine

    engine_->destroyCameraComponent(cameraEntity_);
    utils::EntityManager::get().destroy(cameraEntity_);
    
    if (swapChain_) engine_->destroy(swapChain_);
    if (mainView_) engine_->destroy(mainView_);
    if (renderer_) engine_->destroy(renderer_);
    
    filament::Engine::destroy(engine_);
}

bool MasterRenderer::initialize(void* nativeWindowHandle, int width, int height) {
    width_ = width; 
    height_ = height; 
    nativeWindowHandle_ = nativeWindowHandle;

    engine_ = filament::Engine::create(filament::Engine::Backend::VULKAN);
    if (!engine_) return false;

    renderer_ = engine_->createRenderer();
    mainClearOptions_.clear = true; 
    mainClearOptions_.discard = false; 
    mainClearOptions_.clearColor = {0.03, 0.04, 0.07, 1.0};
    renderer_->setClearOptions(mainClearOptions_);

    swapChain_ = engine_->createSwapChain(nativeWindowHandle_);
    mainView_ = engine_->createView();
    
    cameraEntity_ = utils::EntityManager::get().create();
    camera_ = engine_->createCamera(cameraEntity_);
    
    mainView_->setCamera(camera_);
    mainView_->setAntiAliasing(filament::View::AntiAliasing::FXAA);
    mainView_->setPostProcessingEnabled(true);

    filament::View::BloomOptions bloom = mainView_->getBloomOptions();
    bloom.enabled = true;
    mainView_->setBloomOptions(bloom);

    filament::View::AmbientOcclusionOptions ao;
    ao.enabled = true;
    ao.intensity = 2.0f;      
    ao.radius = 0.3f;          
    ao.bias = 0.005f;         
    ao.power = 1.0f;
    ao.quality = filament::View::QualityLevel::HIGH; 
    
    mainView_->setAmbientOcclusionOptions(ao);
    
    // Color Grading
    // auto toneMapper = new filament::PBRNeutralToneMapper();
    // colorGrading_ = filament::ColorGrading::Builder()
    //     .toneMapper(toneMapper)
    //     .quality(filament::ColorGrading::QualityLevel::HIGH)
    //     .build(*engine_);
    // if (colorGrading_) mainView_->setColorGrading(colorGrading_);
    // delete toneMapper;

    // Input/Camera Controllers
    cameraController_ = std::make_unique<camera::CameraController>(inputState_);
    
    entityPicker_ = std::make_unique<EntityPicker>(*engine_, width, height, camera_, inputState_);

    resize(width, height);
    return true;
}

void MasterRenderer::resize(int width, int height) {
    width_ = width; 
    height_ = height;
    filament::Viewport vp{0u, 0u, static_cast<uint32_t>(width), static_cast<uint32_t>(height)};
    mainView_->setViewport(vp);
    
    const double aspect = static_cast<double>(width) / static_cast<double>(height);
    camera_->setProjection(60.0, aspect, 0.001, 100.0, filament::Camera::Fov::VERTICAL);
    
    if (entityPicker_) {
        entityPicker_->resize(width, height, camera_);
    }
}

bool MasterRenderer::render(GameWorld* world, double elapsedSeconds) {
    if (!renderer_ || !world) return false;

    // 1. Update Camera
    if (cameraController_) {
        cameraController_->update(elapsedSeconds);
        cameraController_->applyTo(camera_);
    }

    // 2. Begin Frame
    if (!renderer_->beginFrame(swapChain_)) {
        return false; 
    }

    // 3. Update Entity Picker and Material Highlighting
    if (world->areMaterialsReady() && entityPicker_) {
        entityPicker_->setScene(world->getPickingScene());
        entityPicker_->update(static_cast<uint32_t>(inputState_->mouseX), static_cast<uint32_t>(inputState_->mouseY),  *renderer_);        

        uint32_t currentHoveredId = entityPicker_->getGrabState()->hoveredEntityId;

        auto& rm = engine_->getRenderableManager();

//         if (lastHoveredEntityId != 0 && lastHoveredEntityId != currentHoveredId) {
//             utils::Entity oldEntity = world->getRenderableEntityForPickId(lastHoveredEntityId);
//             auto rInstance = rm.getInstance(oldEntity);
//             if (rInstance) {
//                 rm.setMaterialInstanceAt(rInstance, 0, world->getBaseMaterialInstance());
//             }
//         }
// 
//         if (currentHoveredId != 0 && currentHoveredId != lastHoveredEntityId) {
//             utils::Entity newEntity = world->getRenderableEntityForPickId(currentHoveredId);
//             auto rInstance = rm.getInstance(newEntity);
//             if (rInstance) {
//                 rm.setMaterialInstanceAt(rInstance, 0, world->getHighlightMaterialInstance());
//             }
//         }

        lastHoveredEntityId = currentHoveredId;
    }

    // 4. Render Main Scene
    if (world->areMaterialsReady()) {
        mainView_->setScene(world->getMainScene());
        renderer_->setClearOptions(mainClearOptions_);
        renderer_->render(mainView_);
    }

    renderer_->endFrame();
    return true;
}

} // namespace filament_example::render