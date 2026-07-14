#pragma once

#include "../SceneBase.hpp"
#include "misc/spline.hpp"
#include <vector>

using namespace cardillo;

class NewtonsCradleScene : public SceneBase {
public:
    const char* sceneName() const override { return "newtons_cradle"; }

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        engine.setGravity(Vector3r(0.0, 0.0, -9.81));

        const real_t t = 0.02;         // Half-thickness for pillars and bars
        const real_t frameHx = 0.4;    // Half-length of the cradle (X axis)
        const real_t frameHy = 0.2;    // Half-width of the cradle (Y axis)
        const real_t frameHz = 0.3;    // Half-height of the pillars (Z axis)

        const real_t topZ = 2.0 * frameHz; // Absolute top Z coordinate
        const real_t bottomZ = 0.0;        // Absolute bottom Z coordinate

        const Vector3r pillarExtents(t, t, frameHz);
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r(-frameHx, -frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r( frameHx, -frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r(-frameHx,  frameHy, frameHz)));
        engine.addStaticBody(physics::CubeShape(pillarExtents), physics::RigidState(Vector3r( frameHx,  frameHy, frameHz)));

        const Vector3r horizontalExtents(frameHx + t, t, t);
        const Vector3r leftPos = Vector3r(0.0, -frameHy, topZ);
        const Vector3r rightPos = Vector3r(0.0,  frameHy, topZ);
        auto leftTopBar = engine.addStaticBody(physics::CubeShape(horizontalExtents), physics::RigidState(leftPos));
        auto rightTopBar = engine.addStaticBody(physics::CubeShape(horizontalExtents), physics::RigidState(rightPos));

        engine.addStaticBody(physics::CubeShape(Vector3r(t, frameHy, t)), physics::RigidState(Vector3r(-frameHx, 0.0, bottomZ + t)));
        engine.addStaticBody(physics::CubeShape(Vector3r(t, frameHy, t)), physics::RigidState(Vector3r( frameHx, 0.0, bottomZ + t)));

        const int numSpheres = 7;
        const real_t radius = 0.038;
        const real_t spacing = 2.01 * radius; // Ensure they are exactly touching
        const real_t ballZ = 0.15;            // Hanging height, safely above the floor

        physics::RigidProps sphereProps = physics::RigidProps::withDensity(7800.0); // Steel density
        sphereProps.restitution_normal = 0.999;
        sphereProps.restitution_tangential = 0.1; 

        for (int i = 0; i < numSpheres; ++i) {
            real_t xPos = (i - (numSpheres-1) / 2.0) * spacing;
            
            physics::RigidState state;
            state.position = Vector3r(xPos, 0.0, ballZ);

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

            Vector3r leftAnchor(xPos, -frameHy, topZ);
            Vector3r rightAnchor(xPos, frameHy, topZ);

            // LINEAR DISTANCE CONSTRAINTS 

            // engine.addLinearDistanceConstraint(sphereEnt, leftTopBar, Vector3r::Zero(), leftAnchor - leftPos);  // Left wire
            // engine.addLinearDistanceConstraint(sphereEnt, rightTopBar, Vector3r::Zero(), rightAnchor - rightPos); // Right wire

            // WIRES BETWEEN SPHERE AND TOP BARS

            auto rightSpline = misc::LinearSpline({rightAnchor, state.position});
            auto leftSpline = misc::LinearSpline({leftAnchor, state.position});
            auto crossSection = physics::BeamCrossSection(0.002, 0.002);
            auto springs = physics::BeamSpringParams::fromMaterial(1e9, 0.3, 500.0);

            auto rightBeam = engine.createBeam(
                rightSpline,
                crossSection,
                springs,
                physics::RigidState(Vector3r::Zero(), Quaternion4r::Identity()),
                physics::RigidProps::withDensity(50.0),
                10
            );

            engine.disableCollisionBetween(sphereEnt, rightBeam.second);
            engine.addTranslationalConstraint(sphereEnt, rightBeam.second, physics::JointFrame{rightBeam.second});
            engine.disableCollisionBetween(rightTopBar, rightBeam.first);
            engine.addTranslationalConstraint(rightTopBar, rightBeam.first, physics::JointFrame{rightBeam.first});

            auto leftBeam = engine.createBeam(
                leftSpline,
                crossSection,
                springs,
                physics::RigidState(Vector3r::Zero(), Quaternion4r::Identity()),
                physics::RigidProps::withDensity(50.0),
                10
            );

            engine.disableCollisionBetween(sphereEnt, leftBeam.second);
            engine.addTranslationalConstraint(sphereEnt, leftBeam.second, physics::JointFrame{ leftBeam.second});
            engine.disableCollisionBetween(leftTopBar, leftBeam.first);
            engine.addTranslationalConstraint(leftTopBar, leftBeam.first, physics::JointFrame{leftBeam.first});
        }
    }
};