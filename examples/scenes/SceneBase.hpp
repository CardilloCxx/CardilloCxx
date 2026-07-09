#pragma once

#include "../src/physics/api/physics.hpp"
#include "../src/misc/types.hpp"
#include "../src/config/path.hpp"

// Abstract base for example scenes. Derive and implement `populate` which
// receives an already-constructed PhysicsEngine to add obstacles and bodies to.
class SceneBase {
public:
    SceneBase() = default;
    virtual ~SceneBase() = default;

    // Short identifier used for config selection and output naming
    virtual const char* sceneName() const = 0;

    // Populate the provided physics engine (add obstacles, bodies, etc.).
    virtual void populate(cardillo::physics::PhysicsEngine& engine) { (void)engine; }

    // Optional per-frame scene update. Called once each simulation step with
    // current time t and timestep dt. Default is no-op.
    virtual void updateScene(cardillo::physics::PhysicsEngine& engine, real_t /*t*/, real_t /*dt*/) { (void)engine; }
};
