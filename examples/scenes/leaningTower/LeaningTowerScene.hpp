#pragma once

#include <iostream>
#include "../SceneBase.hpp"

// Leaning tower scene: stack of slender blocks forming a tower with a small
// lateral offset. Two configurations are used in practice: one with a
// conservative offset where the projection of the centre of mass stays
// safely inside the support polygon, and one where the tower is tilted
// closer to the edge to test stability and collapse.
class LeaningTowerScene : public SceneBase {
public:
    const char* sceneName() const override { return "leaningTower"; }
    LeaningTowerScene() = default;
    ~LeaningTowerScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

        // Slender tower block half extents (x,y horizontal, z vertical)
        const Vector3r half(0.1, 0.03, 0.01);
        const real_t density = (real_t)800.0;
        const int nLevels = 24;

        // Horizontal offset per level. The actual magnitude is controlled
        // via the configuration ("safe" vs "near-edge") by choosing a
        // different epsilon in the scene config files.
        const real_t eps = (real_t)0.00001; // fine-tuned per config

        PhysicsSystem::CubeShape blockShape{half};

        real_t x = 0;

        for (int k = 0; k < nLevels; ++k) {
            const real_t z = half.z() + (real_t)1.999 * half.z() * (real_t)k;
            x += (half.x() / (nLevels - k) ) - eps;

            auto state = PhysicsSystem::RigidState( Vector3r(x, 0.0, z) );
            auto props = PhysicsSystem::RigidProps::withDensity(density);

            auto entity = sys.addRigidBody(blockShape, state, props);
            if (k == 0) sys.makeStatic(entity);
        }

        x = 0.5;
        for (int k = 0; k < nLevels; ++k) {
            const real_t z = half.z() + (real_t)1.999 * half.z() * (real_t)k;
            x += (half.x() / (nLevels - k) ) + eps;

            auto state = PhysicsSystem::RigidState( Vector3r(x, 0.0, z) );
            auto props = PhysicsSystem::RigidProps::withDensity(density);

            auto entity = sys.addRigidBody(blockShape, state, props);
            if (k == 0) sys.makeStatic(entity);
        }
    }
};
