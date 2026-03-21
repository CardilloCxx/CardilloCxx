#pragma once

#include <entt/entt.hpp>

#include "world.hpp"

namespace cardillo {
namespace physics {

class BeamElementUpdater {
public:
    static void updateBeamElementEntity(World& world, entt::entity e);
    static void updateEntities(World& world);
};

} // namespace physics
} // namespace cardillo
