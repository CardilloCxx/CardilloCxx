#pragma once

#include "../src/physics/physics_system.hpp"
#include "../src/misc/types.hpp"

// Abstract base for example scenes. Derive and implement `populate` which
// receives an already-constructed PhysicsSystem to add obstacles and bodies to.
class SceneBase {
public:
    SceneBase() = default;
    virtual ~SceneBase() = default;

    // Populate the provided physics system (add obstacles, bodies, etc.).
    virtual void populate(cardillo::PhysicsSystem& sys) = 0;

    // Optional per-frame scene update. Called once each simulation step with
    // current time t and timestep dt. Default is no-op.
    virtual void updateScene(cardillo::PhysicsSystem& sys, real_t /*t*/, real_t /*dt*/) { (void)sys; }
};
