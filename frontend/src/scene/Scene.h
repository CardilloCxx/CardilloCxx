#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>

#include <filament/Engine.h>
#include <filament/IndexBuffer.h>
#include <filament/RenderableManager.h>
#include <filament/Scene.h>
#include <filament/MaterialInstance.h>
#include <filament/VertexBuffer.h>
#include <utils/Entity.h>
#include <utils/EntityManager.h>
#include <math/mat4.h>

#include "world/GameWorld.h"

namespace frontend::scene {
    void setupCornellBox(GameWorld& world);
} // namespace filament_example::scene
