#pragma once

#include "../SceneBase.hpp"

using namespace cardillo;

class TableScene : public SceneBase {
   public:
    const char* sceneName() const override { return "table"; }

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        engine.setGravity(Vector3r(0, 0, -9.81));

        // Single static cube used as the table.
        const Vector3r tableHalfExtents((real_t)0.45, (real_t)0.30, (real_t)0.035);
        const Vector3r tableCenter((real_t)0.0, (real_t)0.0, (real_t)0.70);
        engine.addStaticBody(physics::CubeShape(tableHalfExtents), physics::RigidState(tableCenter));

        // Thin book-like rigid body, broad face horizontal (flat).
        const Vector3r bookHalfExtents((real_t)0.14, (real_t)0.10, (real_t)0.012);
        const real_t tableTopZ = tableCenter.z() + tableHalfExtents.z();
        const Vector3r bookStart((real_t)0.10, (real_t)0.0, tableTopZ + (real_t)0.045);

        physics::RigidState state;
        state.position = bookStart;
        state.linearVelocity = Vector3r((real_t)1.85, (real_t)0.0, (real_t)-0.05);
        state.angularVelocity = Vector3r((real_t)0.0, (real_t)-1.2, (real_t)0.3);

        physics::RigidProps props((real_t)1.2);
        props.friction = (real_t)0.25;

        engine.addRigidBody(physics::CubeShape(bookHalfExtents), state, props);
    }
};
