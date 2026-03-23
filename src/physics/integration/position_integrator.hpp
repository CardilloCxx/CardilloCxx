#pragma once

#include "../synchronization/derived_entity_sync.hpp"

namespace cardillo {
namespace physics {

class PositionIntegrator {
public:
    static void explicitPositionUpdate(World& world, real_t h);
    static void linearImplicitPositionUpdate(World& world, real_t h);
};

} // namespace physics
} // namespace cardillo
