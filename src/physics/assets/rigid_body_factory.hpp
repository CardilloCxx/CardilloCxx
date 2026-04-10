#pragma once

#include <entt/entt.hpp>

#include "../world.hpp"

namespace cardillo::physics {

struct RigidBodyFactory {
    static entt::entity create(World& system, const physics::RigidShape& shape, const physics::RigidState& state, const physics::RigidProps& props);
};

}  // namespace cardillo::physics
