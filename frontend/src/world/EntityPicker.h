#pragma once

#include <cstdint>
#include <memory>

#include <filament/Camera.h>
#include <filament/Engine.h>
#include <filament/RenderTarget.h>
#include <filament/Renderer.h>
#include <filament/Scene.h>
#include <filament/Texture.h>
#include <filament/View.h>
#include <math/vec4.h>

#include "platform/Platform.h"

namespace frontend {

class GameWorld; 

class EntityPicker {

public:
    struct GrabState{
        bool isGrabbing = false;
        uint32_t hoveredEntityId = 0;
        uint32_t hoveredEntityRenderableId = 0;

        filament::math::float3 inBodyGrabPoint;
        filament::math::float3 inWorldGrabPoint;

        float initialGrabDistance = 0.0f;
        filament::math::float3 mouseWorldPos;
    };

    struct Impl;

    EntityPicker(filament::Engine& engine, int width, int height, filament::Camera* camera, const platform::InputState* inputState);
    ~EntityPicker();

    void resize(int width, int height, filament::Camera* camera);
    void setCamera(filament::Camera* camera);
    void setScene(filament::Scene* scene);
    void setGameWorld(frontend::GameWorld* gameWorld);

    uint32_t getHoveredEntityId() const;
    const GrabState* getGrabState() const { return grabState_.get(); };
    void update(uint32_t x, uint32_t y, filament::Renderer& renderer);
    void updateGrabbed(uint32_t x, uint32_t y, filament::math::float3 worldPos, uint32_t entityId, filament::Renderer& renderer, filament::Camera& camera);

    filament::Scene* scene() const noexcept;
    filament::View* view() const noexcept;
    filament::RenderTarget* renderTarget() const noexcept;

private:
    const platform::InputState* input_ = nullptr;
    frontend::GameWorld* gameWorld_ = nullptr;

    std::unique_ptr<Impl> impl_;
    std::unique_ptr<GrabState> grabState_;
};

filament::math::float4 encodeEntityPickColor(uint32_t entityId);

} // namespace filament_example::render
