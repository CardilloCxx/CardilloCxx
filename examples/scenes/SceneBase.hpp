#pragma once

#include "../src/physics/physics_system.hpp"

// Abstract base for example scenes. Derive and implement `populate` which
// receives an already-constructed PhysicsSystem to add obstacles and bodies to.
class SceneBase {
public:
    SceneBase() = default;
    virtual ~SceneBase() = default;

    // Populate the provided physics system (add obstacles, bodies, etc.).
    virtual void populate(cardillo::PhysicsSystem& sys) = 0;
};
