#pragma once

#include "../world.hpp"

namespace cardillo {
namespace physics {

class Trajectory {
   public:
    // Sets the pose of every kinematically driven body to its prescribed value at
    // t_n + evalOffset (t_n = elapsed time at the beginning of the step) and its velocity to the
    // mean velocity over the step, (pose(t_n + dt) - pose(t_n)) / dt (angular velocity in the body
    // frame). The Moreau-theta integrator uses evalOffset = (1 - theta) dt, i.e. the
    // time of the intermediate configuration q_{n+theta}.
    static void update(World& world, real_t dt, real_t evalOffset = 0);
};

}  // namespace physics
}  // namespace cardillo
