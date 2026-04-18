#pragma once

#include "../world.hpp"

namespace cardillo {
namespace physics {

class Trajectory {
   public:
    static void update(World& world, real_t dt);
};

}  // namespace physics
}  // namespace cardillo
