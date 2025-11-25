#pragma once

#include "../SceneBase.hpp"
#include <vector>

using namespace cardillo;

class SoftbodyTestScene : public SceneBase {
public:
    const char* sceneName() const override { return "softbody"; }
    SoftbodyTestScene() = default;
    ~SoftbodyTestScene() override = default;

    void populate(cardillo::PhysicsSystem& sys) override {
        using namespace cardillo;

    PhysicsSystem::PlaneShape groundShape{Vector3r(0,0,1), Vector3r(0,1,0), (real_t)15.0, (real_t)15.0};
    PhysicsSystem::RigidState groundState; groundState.position = Vector3r(0.0,0.0,0.0); groundState.orientation = Quaternion4r::Identity();
    (void)sys.addStaticBody(groundShape, groundState);

        // Load soft body from OBJ and place it 0.7 m above the floor (floor top at z = 0)
        const std::string path = "res/meshes/teapot.obj";
        const real_t k = (real_t)200000.0;
        const real_t d = (real_t)3000.0;
        const Vector3r position(0.0,0.0,0.01);
        const Quaternion4r orientation = Quaternion4r(Eigen::AngleAxis<real_t>(M_PI_2, Vector3r::UnitX()));
        const Vector3r v = Vector3r::Zero();
        const Vector3r w = Vector3r::Zero();
        const real_t totalMass = (real_t)2.0; // distribute 2 kg uniformly
        (void)sys.addSoftBody(path, k, d, position, orientation, v, w, totalMass);

        // Add a cube to crush the soft body
        {
            PhysicsSystem::CubeShape shape{Vector3r(0.2,0.2,0.2)}; PhysicsSystem::RigidState st; st.position = Vector3r(0.0,0.0,2.5); st.orientation = Quaternion4r::Identity(); PhysicsSystem::RigidProps pr; pr.mass = (real_t)100.0; sys.addRigidBody(shape, st, pr);
        }
    }
};
