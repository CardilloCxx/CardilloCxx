#pragma once

#include "../SceneBase.hpp"
#include <vector>

using namespace cardillo;

class NewtonsCradleScene : public SceneBase {
public:
    const char* sceneName() const override { return "newtons_cradle"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        // Set standard gravity, matching standard Z-up orientation
        engine.setGravity(Vector3r(0.0, 0.0, -9.81));

        // --- 1. Frame Dimensions (Half-Extents) ---
        const real_t t = 0.02;         // Half-thickness for pillars and bars
        const real_t frameHx = 0.4;    // Half-length of the cradle (X axis)
        const real_t frameHy = 0.2;    // Half-width of the cradle (Y axis)
        const real_t frameHz = 0.3;    // Half-height of the pillars (Z axis)

        const real_t topZ = 2.0 * frameHz; // Absolute top Z coordinate
        const real_t bottomZ = 0.0;        // Absolute bottom Z coordinate

        // --- 2. Build the Frame (8 Cubes) ---
        // 4 Pillars
        const Vector3r pillarExtents(t, t, frameHz);
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r(-frameHx, -frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r( frameHx, -frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r(-frameHx,  frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r( frameHx,  frameHy, frameHz)));

        // 2 Top bars going across (connecting pillars along the length/X-axis)
        // We add '+ t' to the X extent so it caps perfectly over the pillars
        const Vector3r horizontalExtents(frameHx + t, t, t);
        const Vector3r leftPos = Vector3r(0.0, -frameHy, topZ);
        const Vector3r rightPos = Vector3r(0.0,  frameHy, topZ);
        auto leftTopBar = engine.addStaticBody(physics::CubeShape(horizontalExtents), physics::RigidState(leftPos));
        auto rightTopBar = engine.addStaticBody(physics::CubeShape(horizontalExtents), physics::RigidState(rightPos));

        // 2 Bottom foots connecting the bottoms
        engine.addStaticBody(physics::CubeShape(Vector3r(t, frameHy, t)), physics::RigidState(Vector3r(-frameHx, 0.0, bottomZ + t)));
        engine.addStaticBody(physics::CubeShape(Vector3r(t, frameHy, t)), physics::RigidState(Vector3r( frameHx, 0.0, bottomZ + t)));


        // --- 3. Spheres Setup ---
        const int numSpheres = 5;
        const real_t radius = 0.038;
        const real_t spacing = 2.01 * radius; // Ensure they are exactly touching
        const real_t ballZ = 0.15;           // Hanging height, safely above the floor

        // Elastic properties for ideal energy transfer
        physics::RigidProps sphereProps = physics::RigidProps::withDensity(7800.0); // Steel density
        sphereProps.restitution_normal = 0.999;
        sphereProps.restitution_tangential = 0.1; 

        for (int i = 0; i < numSpheres; ++i) {
            // Center the 5 spheres around the origin (x = 0)
            real_t xPos = (i - 2) * spacing;
            
            physics::RigidState state;
            state.position = Vector3r(xPos, 0.0, ballZ);

            // Give the outer sphere (leftmost) initial velocity to the side
            if (i == 0) {
                state.linearVelocity = Vector3r(-2.0, 0.0, 0.0);
            }
            if (i == 1) {
                state.linearVelocity = Vector3r(-2.0, 0.0, 0.0);
            }
            if (i == numSpheres - 1) {
                state.linearVelocity = Vector3r(2.0, 0.0, 0.0);
            }

            auto sphereEnt = engine.addRigidBody(physics::SphereShape(radius), state, sphereProps);

            // --- 4. V-Shaped Wires ---
            // A true Newton's Cradle uses rigid distance constraints to prevent twisting.
            // Assuming your engine supports standard distance constraints:
            Vector3r leftAnchor(xPos, -frameHy, topZ);
            Vector3r rightAnchor(xPos, frameHy, topZ);

            engine.addLinearDistanceConstraint(sphereEnt, leftTopBar, Vector3r::Zero(), leftAnchor - leftPos);  // Left wire
            engine.addLinearDistanceConstraint(sphereEnt, rightTopBar, Vector3r::Zero(), rightAnchor - rightPos); // Right wire
        }
    }
};