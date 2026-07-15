#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Domino scene: builds a layered domino tower near the origin.
class DominoScene : public SceneBase {
public:
    const char* sceneName() const override { return "domino"; }
    DominoScene() = default;
    ~DominoScene() override = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        // Create a ground plate (static cube via unified API)
        physics::CubeShape groundShape{Vector3r(3.0, 3.0, 0.05)};
        physics::RigidState groundState;
        groundState.position = Vector3r(0.0, 0.0, -0.05);

        // ground as free rigid body with constraints
        // auto ground = engine.addRigidBody(groundShape, groundState, physics::RigidProps::withDensity(100.0));
        // engine.addTranslationalConstraint(ground, entt::null, physics::JointFrame(ground));
        // engine.addRotationConstraint(ground, entt::null, physics::JointFrame(ground));

        // alternative: ground as static body
        engine.addStaticBody(groundShape, groundState);

        // Domino dims: x=length/2, y=thickness/2, z=height/2
        const Vector3r dominoHalf((real_t)0.024, (real_t)0.00375, (real_t)0.012); // length 9.6cm, thickness 1.5cm, height 4.8cm
        const real_t density = (real_t)800.0;
        const int layers = 11;  // 186; //11; // 26; //53;
        const int gridN = 8;    // 8; //16; // 16; //32;
        const Vector3r baseCenter(0.0, 0.0, 0.0);
        const real_t gapLong = (real_t)0.004; // small longitudinal spacing
        const real_t extraLayerGap = (real_t)-0.0001;

        spawnDominoTowerStructure(engine, layers, gridN, dominoHalf, density, baseCenter, gapLong, extraLayerGap, true);

        // sphere falling on the tower
        // auto sphere = engine.addRigidBody(
        //     physics::SphereShape{0.05}, 
        //     physics::RigidState{
        //         Vector3r(0.0, 0.0, 1.5), // position
        //         Vector3r(0.0, 0.0, -1.0), // linear velocity
        //         Vector3r(5.0, 0.0, 0.0)}, // angular velocity
        //     physics::RigidProps::withDensity(19250.0)
        // );

        // sphere crashing the tower
        // auto sphere = engine.addRigidBody(
        //     physics::SphereShape{0.05}, 
        //     physics::RigidState{
        //         Vector3r(0.3, 0.0, 0.1), // position
        //         Vector3r(-5.0, 0.0, 0.0), // linear velocity
        //         Vector3r(0.0, 0.0, 0.0)}, // angular velocity
        //     physics::RigidProps::withDensity(19250.0)
        // );

        // auto helixSpline = misc::HelixSpline(Vector3r(0.0, 0.0, 1.3), -Vector3r::UnitZ(), (real_t)0.1, (real_t)0.15, (real_t)10.0, Vector3r::UnitX());

        // engine.addTrajectory(
        //     sphere,
        //     [helixSpline](real_t t) {
        //         TrajectoryPose pose;
        //         pose.first = helixSpline.sample(std::fmod(t, (real_t)1.0)).position;
        //         return pose;
        //     },
        //     std::nullopt);
    }

private:
    // local mass helper (same logic as examples/rigid_body_test)
    static inline real_t massFromDensity(const Vector3r& halfExtents, real_t density, real_t minMass = (real_t)0.05) {
        const real_t volume = (real_t)8.0 * halfExtents.x() * halfExtents.y() * halfExtents.z();
        return std::max(minMass, density * volume);
    }

    void spawnDominoTowerStructure(
        physics::PhysicsEngine& engine,
        int layers,
        int N,
        const Vector3r& half,
        real_t density,
        const Vector3r& baseCenter,
        real_t gapLong = (real_t)0.002,
        real_t extraLayerGap = (real_t)0.0,
        bool initialImpulse = true
    ) {
        (void)gapLong; (void)extraLayerGap; // parameters kept for compatibility with older callers
        if (layers <= 0) return;
        const real_t L = (real_t)2.0 * half.x(); // long edge (in-plane)
        const real_t W = (real_t)2.0 * half.y(); // thickness (in-plane)
        const real_t H = (real_t)2.0 * half.z(); // height (upright)
        const real_t Offset = (real_t) 0.3 * W;

        // Grid spacing is domino length minus thickness, per specification
        const real_t s = std::max<real_t>((real_t)1e-6, L - (W - Offset));
        const int Ncells = std::max(2, (N / 2) * 2); // ensure even

        // A lambda function that places a domino between (i,j,k) and (i, j+1, k) or (i+1, j, k)
        auto placeDomino = [&](int i, int j, int k, bool alongY) {
            Vector3r c = baseCenter;
            // compute position offsets
            const real_t offsetX = ((real_t)i - (real_t)(Ncells - 1) * (real_t)0.5) * s - (0.5 * Offset) * alongY;
            const real_t offsetY = ((real_t)j - (real_t)(Ncells - 1) * (real_t)0.5) * s - (0.5 * Offset) * !alongY;
            c.x() += offsetX;
            c.y() += offsetY;
            const real_t z = baseCenter.z() + half.z() + (real_t)k * ( (real_t)2.0 * half.z() + extraLayerGap );
            c.z() = z;  
            physics::CubeShape shape{half};
            Quaternion4r q = Quaternion4r::Identity();
            if (alongY) {
                // rotate 90 deg about Z to align long axis along Y
                q = Quaternion4r(Eigen::AngleAxis<real_t>((real_t)M_PI_2, Vector3r::UnitZ()));
                c.x() -= (L/2.0 - W/2.0);
                c.y() += (L/2.0 - W/2.0);
            }
            const real_t m = massFromDensity(half, density);

            Vector3r vel = Vector3r::Zero();
            if (initialImpulse && (i == Ncells -1 && (j == Ncells /2 || j == Ncells /2 - 1) && k == layers -4)) {
                vel = Vector3r(4.0, 0.0, -1.0) * 2;
            }

            physics::RigidState state; state.position = c; state.orientation = q; state.linearVelocity = vel; state.angularVelocity = Vector3r::Zero();
            physics::RigidProps props; props.mass = m;
            engine.addRigidBody(shape, state, props);

        };

        // Place layers 
        for (int layer = 0; layer < layers; ++layer) {
            const int off = 0;
            const int k = layer;
            // Each layer is a grid of dominos
            if((layer + off) % 4 == 0) 
            {
                // place all parrallel to x from (i, j) to (i, j+1) and (i, j+2) to (i, j+3) and so on
                for (int i = 0; i < Ncells; ++i) {
                    for (int j = 0; j < Ncells; j +=2) {
                        placeDomino(i, j, k, true);
                    }
                }
            }else if ((layer + off) % 4 == 1)
            {
                // First place parralel along y from (i, j) to (i+1, j) and (i+2, j) to (i+3, j) and so on at bottom and top only:
                for(int i = 0; i < Ncells; i +=2)
                {
                    placeDomino(i, 0, k, false);
                    placeDomino(i, Ncells -1, k, false);
                }
                // Then place parralel along x along at left and right only:
                for(int j = 1; j < Ncells -1; j +=2)
                {
                    placeDomino(0, j, k, true);
                    placeDomino(Ncells -1, j, k, true);
                }
                // now we fill the middle parralel to y from (i+1, j+1) to (i+2, j+1) and so on
                for(int i =1; i < Ncells -1; i +=2)
                {
                    for(int j =1; j < Ncells -1; j +=1)
                    {
                        placeDomino(i, j, k, false);
                    }
                }
            }else if ((layer + off) % 4 == 2)
            {
                // place parallel to x left and right only:
                for(int j = 0; j < Ncells; j +=2)
                {
                    placeDomino(0, j, k, true);
                    placeDomino(Ncells -1, j, k, true);
                }
                // place parallel to y top and bottom only:
                for(int i = 1; i < Ncells -1; i +=2)
                {
                    placeDomino(i, 0, k, false);
                    placeDomino(i, Ncells -1, k, false);
                }
                // now fill the middle
                for(int i =1; i < Ncells -1; i +=1)
                {
                    for(int j =1; j < Ncells -1; j +=2)
                    {
                        placeDomino(i, j, k, true);
                    }
                }
            }else if ((layer + off) % 4 == 3)
            {
                for(int i = 0; i < Ncells; i +=2)
                {
                for(int j = 0; j < Ncells; j ++)
                {
                        placeDomino(i, j, k, false);
                }
                }
            }
        }
    }
};
