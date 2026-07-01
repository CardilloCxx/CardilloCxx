#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

// Jenga scene: builds a layered block tower using the existing helper logic.
class JengaScene : public SceneBase {
public:
    const char* sceneName() const override { return "jenga"; }
    JengaScene() = default;
    ~JengaScene() override = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

    // Ground (static cube via unified API)
    physics::CubeShape groundShape{Vector3r(15.0, 15.0, 0.5)};
    physics::RigidState groundState; groundState.position = Vector3r(0.0, 0.0, -0.5); groundState.orientation = Quaternion4r::Identity();
    engine.addStaticBody(groundShape, groundState);

        // Build a Jenga tower
        const Vector3r blockHalf((real_t)0.075, (real_t)0.0225, (real_t)0.0125); // example block half extents
        const int layers = 12;
        const real_t gap = (real_t)0.002;
        const real_t density = (real_t)600.0;
        const Vector3r baseCenter(0.0, 0.0, 0.0);
        const int blocksPerLayer = 3;
        const real_t extraLayerGap = (real_t)0.0;

        // place tower using same algorithm as previous helper
        const Quaternion4r q0 = Quaternion4r::Identity();
        const Quaternion4r q90 = Quaternion4r(Eigen::AngleAxis<real_t>((real_t)M_PI_2, Vector3r::UnitZ()));
        const real_t fullY = (real_t)2.0 * blockHalf.y();
        const real_t fullZ = (real_t)2.0 * blockHalf.z();
        const real_t baseZ = baseCenter.z() + blockHalf.z();

        for (int layer = 0; layer < layers; ++layer) {
            const bool alongX = (layer % 2 == 0);
            const Quaternion4r q = alongX ? q0 : q90;
            const real_t rowWidth = blocksPerLayer * fullY + (blocksPerLayer - 1) * gap;
            const real_t firstOffset = -rowWidth * (real_t)0.5 + fullY * (real_t)0.5;
            const real_t z = baseZ + (real_t)layer * (fullZ + extraLayerGap);
            for (int i = 0; i < blocksPerLayer; ++i) {
                Vector3r c = baseCenter;
                const real_t step = fullY + gap;
                const real_t offset = firstOffset + (real_t)i * step;
                if (alongX) c = Vector3r(baseCenter.x(), baseCenter.y() + offset, z);
                else c = Vector3r(baseCenter.x() + offset, baseCenter.y(), z);
                physics::CubeShape blkShape{blockHalf};
                physics::RigidState blkState; blkState.position = c; blkState.orientation = q;
                physics::RigidProps blkProps; blkProps.mass = std::max((real_t)0.05, density * (real_t)8.0 * blockHalf.x() * blockHalf.y() * blockHalf.z());
                engine.addRigidBody(blkShape, blkState, blkProps);
            }
        }
    }
};
