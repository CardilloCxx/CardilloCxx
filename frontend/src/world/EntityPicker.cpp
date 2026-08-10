#include "world/EntityPicker.h"
#include "world/GameWorld.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>

#include <filament/Renderer.h>
#include <filament/Texture.h>
#include <filament/Viewport.h>

namespace frontend {

struct ReadbackState {
    std::array<float, 4> buffer{};
    bool isWaitingForGPU = false;
    bool isNewDataReady = false;

    filament::math::float3 worldPos;
    uint32_t entityID;
};

void onReadback(void* buffer, size_t size, void* user) {
    auto* state = static_cast<ReadbackState*>(user);
    if (!state) return;
    
    state->isNewDataReady = true;
    state->isWaitingForGPU = false;

    state->worldPos = {static_cast<float*>(buffer)[0], static_cast<float*>(buffer)[1], static_cast<float*>(buffer)[2]};

    float packedId = static_cast<float*>(buffer)[3];
    uint32_t entityID;
    std::memcpy(&entityID, &packedId, sizeof(uint32_t));

    state->entityID = entityID;
}

struct EntityPicker::Impl {
    filament::Engine* engine = nullptr;
    filament::View* view = nullptr;
    filament::Scene* scene = nullptr;
    filament::Texture* colorTexture = nullptr;
    filament::Texture* depthTexture = nullptr;
    filament::RenderTarget* renderTarget = nullptr;
    
    ReadbackState readbackState{};
    uint32_t lastKnownHoveredId = 0u;
    
    int width = 0;
    int height = 0;
};

EntityPicker::EntityPicker(filament::Engine& engine, int width, int height, filament::Camera* camera, const platform::InputState* inputState)
    : impl_(std::make_unique<Impl>()), grabState_(std::make_unique<GrabState>()), input_(inputState) {
    impl_->engine = &engine;
    impl_->width = width;
    impl_->height = height;
    resize(width, height, camera);
}

EntityPicker::~EntityPicker() {
    if (!impl_ || !impl_->engine) return;
    if (impl_->renderTarget) impl_->engine->destroy(impl_->renderTarget);
    if (impl_->depthTexture) impl_->engine->destroy(impl_->depthTexture);
    if (impl_->colorTexture) impl_->engine->destroy(impl_->colorTexture);
    if (impl_->view) impl_->engine->destroy(impl_->view);
    if (impl_->scene) impl_->engine->destroy(impl_->scene);
}

void EntityPicker::resize(int width, int height, filament::Camera* camera) {
    impl_->width = width;
    impl_->height = height;
    if (!impl_->engine) return;

    // Reset state on resize to prevent reading stale/out-of-bounds data
    impl_->readbackState = {}; 
    
    if (impl_->renderTarget) impl_->engine->destroy(impl_->renderTarget);
    if (impl_->depthTexture) impl_->engine->destroy(impl_->depthTexture);
    if (impl_->colorTexture) impl_->engine->destroy(impl_->colorTexture);
    if (impl_->view) impl_->engine->destroy(impl_->view);

    if (!impl_->scene) {
        impl_->scene = impl_->engine->createScene();
    }
    
    impl_->view = impl_->engine->createView();
    impl_->view->setScene(impl_->scene);
    impl_->view->setCamera(camera);
    impl_->view->setPostProcessingEnabled(false);
    impl_->view->setAntiAliasing(filament::View::AntiAliasing::NONE);
    impl_->view->setViewport({0u, 0u, static_cast<uint32_t>(width), static_cast<uint32_t>(height)});

    impl_->colorTexture = filament::Texture::Builder()
        .width(static_cast<uint32_t>(width))
        .height(static_cast<uint32_t>(height))
        .levels(1)
        .usage(filament::Texture::Usage::COLOR_ATTACHMENT | filament::Texture::Usage::SAMPLEABLE | filament::Texture::Usage::BLIT_SRC)
        .format(filament::Texture::InternalFormat::RGBA32F)
        .build(*impl_->engine);
        
    impl_->depthTexture = filament::Texture::Builder()
        .width(static_cast<uint32_t>(width))
        .height(static_cast<uint32_t>(height))
        .levels(1)
        .usage(filament::Texture::Usage::DEPTH_ATTACHMENT)
        .format(filament::Texture::InternalFormat::DEPTH32F)
        .build(*impl_->engine);
        
    impl_->renderTarget = filament::RenderTarget::Builder()
        .texture(filament::RenderTarget::AttachmentPoint::COLOR, impl_->colorTexture)
        .texture(filament::RenderTarget::AttachmentPoint::DEPTH, impl_->depthTexture)
        .build(*impl_->engine);
        
    impl_->view->setRenderTarget(impl_->renderTarget);
}

void EntityPicker::setCamera(filament::Camera* camera) {
    if (impl_->view) impl_->view->setCamera(camera);
}

void EntityPicker::setScene(filament::Scene* scene) {
    impl_->scene = scene;
    if (impl_->view) impl_->view->setScene(scene);
}

void EntityPicker::setGameWorld(frontend::GameWorld* gameWorld) {
    gameWorld_ = gameWorld;
}
void EntityPicker::updateGrabbed(uint32_t x, uint32_t y, filament::math::float3 worldPos, uint32_t entityId, filament::Renderer& renderer, filament::Camera& camera) {
    if(!input_ || !grabState_) return;

    bool grabKeyPressed = input_->mouseLeftDown;
    if (grabKeyPressed == grabState_->isGrabbing) return;

    // Stop Grabbing
    if (!grabKeyPressed && grabState_->isGrabbing) {
        grabState_->isGrabbing = false;
        grabState_->hoveredEntityId = 0;
        return;
    }

    if (!grabKeyPressed) return;

    // Start Grabbing
    if (entityId < 1) return; // No entity to grab

    grabState_->isGrabbing = true;
    const filament::math::float3 rayOrigin = camera.getPosition();
    grabState_->initialGrabDistance = length(worldPos - rayOrigin);
    grabState_->inWorldGrabPoint = worldPos;

    utils::Entity entity = gameWorld_->getRenderableEntityForPickId(entityId);
    grabState_->hoveredEntityId = entityId;
    grabState_->hoveredEntityRenderableId = entity.getId();
    
    filament::TransformManager& transformManager = gameWorld_->getEngine().getTransformManager();
    
    auto transformInstance = transformManager.getInstance(entity);
    if (!transformInstance) return; // Safety check in case entity lacks a transform

    const filament::math::mat4f worldTransform = transformManager.getWorldTransform(transformInstance);
    
    const filament::math::float3 entityPosition = worldTransform[3].xyz;
    const filament::math::mat3f entityRotation = worldTransform.upperLeft();

     filament::math::mat3f pureRotation(
        normalize(entityRotation[0]),
        normalize(entityRotation[1]),
        normalize(entityRotation[2])
    );

    grabState_->inBodyGrabPoint = transpose(pureRotation) * (worldPos - entityPosition);
}

void EntityPicker::update(uint32_t x, uint32_t y, filament::Renderer& renderer) {
    if (!impl_->view || !impl_->renderTarget || !impl_->scene) return;
    auto& camera = impl_->view->getCamera();

    if (impl_->readbackState.isNewDataReady) {
                                 
        impl_->lastKnownHoveredId = impl_->readbackState.entityID;
        impl_->readbackState.isNewDataReady = false;

        updateGrabbed(x, y, impl_->readbackState.worldPos, impl_->lastKnownHoveredId, renderer, camera);
    }

    if (grabState_){
        float ndcX = (static_cast<float>(x) / static_cast<float>(impl_->width)) * 2.0f - 1.0f;
        float ndcY = 1.0f - (static_cast<float>(y) / static_cast<float>(impl_->height)) * 2.0f;

        filament::math::mat4 invViewProj = inverse(camera.getProjectionMatrix() * camera.getViewMatrix());
        filament::math::float4 targetWorld = invViewProj * filament::math::float4(ndcX, ndcY, -1.0f, 1.0f);
        filament::math::float3 origin = camera.getPosition();
        filament::math::float3 direction = normalize(targetWorld.xyz / targetWorld.w - origin);

        grabState_->mouseWorldPos = origin + direction * grabState_->initialGrabDistance;
    }

    if (impl_->readbackState.isWaitingForGPU) {
        return; 
    }

    filament::Renderer::ClearOptions pickClearOptions{};
    pickClearOptions.clear = true;
    pickClearOptions.discard = false;
    pickClearOptions.clearColor = {0.0f, 0.0f, 0.0f, 0.0f};
    renderer.setClearOptions(pickClearOptions);
    renderer.render(impl_->view);
    
    if (impl_->engine) {
        impl_->engine->flush(); // Ensure commands are pushed
    }

    const int clampedX = std::clamp(static_cast<int>(x), 0, impl_->width - 1);
    const int clampedY = std::clamp(static_cast<int>(impl_->height - 1 - static_cast<int>(y)), 0, impl_->height - 1);

    impl_->readbackState.isWaitingForGPU = true;
    
    // Entity ID
    renderer.readPixels(
        impl_->renderTarget, 
        static_cast<uint32_t>(clampedX), 
        static_cast<uint32_t>(clampedY), 
        1u, 1u,
        filament::backend::PixelBufferDescriptor(
            impl_->readbackState.buffer.data(), 
            impl_->readbackState.buffer.size() * sizeof(float),
            filament::backend::PixelDataFormat::RGBA, 
            filament::backend::PixelDataType::FLOAT,
            [](void* buffer, size_t size, void* user) {
                onReadback(buffer, size, user);
            }, 
            &impl_->readbackState
        )
    );
}

uint32_t EntityPicker::getHoveredEntityId() const {
    return impl_->lastKnownHoveredId;
}

filament::Scene* EntityPicker::scene() const noexcept { return impl_->scene; }
filament::View* EntityPicker::view() const noexcept { return impl_->view; }
filament::RenderTarget* EntityPicker::renderTarget() const noexcept { return impl_->renderTarget; }

filament::math::float4 encodeEntityPickColor(uint32_t entityId) {
    const uint32_t encoded = entityId + 1u;
    return {
        static_cast<float>((encoded >> 16) & 0xff) / 255.0f,
        static_cast<float>((encoded >> 8) & 0xff) / 255.0f,
        static_cast<float>(encoded & 0xff) / 255.0f,
        1.0f
    };
}

} // namespace filament_example::render