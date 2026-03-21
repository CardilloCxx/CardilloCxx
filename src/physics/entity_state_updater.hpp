#pragma once

#include "beam_element_updater.hpp"

namespace cardillo {
namespace physics {

class EntityStateUpdater {
public:
    static void explicitPositionUpdate(World& world, real_t h);
    static void linearImplicitPositionUpdate(World& world, real_t h);
};

} // namespace physics
} // namespace cardillo
