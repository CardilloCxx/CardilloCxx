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

    // Optional per-step scene update, called once before each simulation step with the time t at
    // the beginning of the step. Default is no-op. The step size is engine.timeStep().
    virtual void updateScene(cardillo::physics::PhysicsEngine& engine, real_t /*t*/) { (void)engine; }

    // Deprecated: dt is redundant (engine.timeStep()). Kept so that scenes overriding the old
    // signature still work; the default forwards to updateScene(engine, t).
    [[deprecated("override updateScene(engine, t) and use engine.timeStep() if the step size is needed")]] virtual void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t,
                                                                                                                                   real_t /*dt*/) {
        updateScene(engine, t);
    }
};
