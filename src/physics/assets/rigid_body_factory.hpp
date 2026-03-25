#pragma once

#include <entt/entt.hpp>

#include "../world.hpp"

namespace cardillo::physics {

struct RigidBodyFactory {
    static entt::entity create(World& system,
                               const World::RigidShape& shape,
                               const World::RigidState& state,
                               const World::RigidProps& props);
};

} // namespace cardillo::physics
