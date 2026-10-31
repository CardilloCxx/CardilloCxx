#pragma once

#include <entt/entt.hpp>

#include "../world.hpp"

namespace cardillo {
namespace physics {

class DerivedEntitySync {
   public:
    static void updateBeamElementEntity(World& world, entt::entity e);
    // evalOffset: see Trajectory::update().
    static void updateEntities(World& world, real_t dt, real_t evalOffset = 0);
};

}  // namespace physics
}  // namespace cardillo
