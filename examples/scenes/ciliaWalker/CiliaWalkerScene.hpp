#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <vector>

using namespace cardillo;

class CiliaWalkerScene : public SceneBase {
public:
    const char* sceneName() const override { return "cilia_walker"; }
    CiliaWalkerScene() = default;
    ~CiliaWalkerScene() = default;

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace misc;

        engine.setGravity(Vector3r(0, 0, -9.81));

        const real_t L  = (real_t)0.0045;   // cilia length
        const real_t w_base = (real_t)0.0015; // base width
        const real_t w_tip  = (real_t)0.0004; // tip width
        const real_t E  = (real_t)1.2e5;    // Young's modulus
        const real_t nu = (real_t)0.45;     
        const real_t rho = (real_t)750.0; 

        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        physics::RigidProps beamProps = physics::RigidProps::withDensity(rho);
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());

        const size_t segments = 8;

        const auto section = physics::BeamCrossSection::procedural([w_base, w_tip, L](real_t s) {
            real_t current_width = w_base + (w_tip - w_base) * s; 
            return physics::BeamCrossSection::rectangular(current_width, w_base);
        });

        const real_t plateRadius = (real_t)0.2; 
        const real_t plateHeight = (real_t)0.01;
        // physics::CylinderShape floorShape(plateRadius, plateHeight);
        physics::CubeShape floorShape(Vector3r(plateRadius, plateRadius, plateHeight));
        physics::RigidProps floorProps;
        floorProps.mass = std::nullopt; 
        floorProps.visual = true;
        entt::entity vibratingPlate = engine.addRigidBody(floorShape, physics::RigidState(Vector3r(0.0, 0.0, -0.005)), floorProps);

         engine.addTrajectory(vibratingPlate, [](real_t t) {
            TrajectoryPose pose;
            const real_t freq = (real_t)40.0;     
            const real_t gamma = (real_t) 2.6;
            const real_t amplitude = 9.81 * gamma / (4.0 * M_PI * M_PI * freq * freq);
            pose.first = Vector3r(0.0, 0.0, -0.005 + amplitude * std::sin(2.0 * M_PI * freq * t));
            pose.second = Quaternion4r::Identity();
            return pose;
        }, std::nullopt);
        // engine.makeStatic(vibratingPlate);

        const real_t arenaRadius = (real_t)0.12; 
        const real_t wallHeight = (real_t)0.03;
        const int numWallSegments = 32;
        const real_t wallSegmentThickness = (real_t)0.008;

        for (int i = 0; i < numWallSegments; ++i) {
            real_t theta1 = (static_cast<real_t>(i) / numWallSegments) * (2.0 * M_PI);
            real_t theta2 = (static_cast<real_t>(i + 1) / numWallSegments) * (2.0 * M_PI);
            real_t segLen = 2.0 * arenaRadius * std::sin(M_PI / numWallSegments);
            real_t midTheta = (theta1 + theta2) * 0.5;

            physics::CubeShape wallSegment(Vector3r(wallSegmentThickness * 0.5, segLen * 0.5, wallHeight * 0.5));
            Vector3r wallPos(arenaRadius * std::cos(midTheta), arenaRadius * std::sin(midTheta), wallHeight * 0.5);
            
            Quaternion4r wallRot(Eigen::AngleAxis<real_t>(midTheta, Vector3r::UnitZ()));
            engine.addStaticBody(wallSegment, physics::RigidState(wallPos, Vector3r::Zero(), wallRot));
        }

        // Walker's body

        const real_t bodySize = (real_t)0.028; 
        const real_t bodyHeight = (real_t)0.01; 
        physics::CubeShape bodyShape(Vector3r(bodySize * 0.5, bodySize * 0.5, bodyHeight * 0.5));
        physics::RigidProps bodyProps = physics::RigidProps::withDensity(1500.0); 
        
        physics::RigidState bodyState(Vector3r(0.0, 0.0, L + bodyHeight * 0.5 + plateHeight * 0.5), Vector3r(0.04, 0.001, 0.0));
        entt::entity walkerBody = engine.addRigidBody(bodyShape, bodyState, bodyProps);
        std::cout << "Walker Mass: " << engine.getMass(walkerBody) << std::endl;

        const int gridSize = 5;
        const real_t spacing = (real_t)0.006; 
        const real_t gridOffset = (gridSize - 1) * spacing * 0.5;

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                real_t xPos = -gridOffset + i * spacing;
                real_t yPos = -gridOffset + j * spacing;

                Vector3r p0(xPos, yPos, -bodyHeight * 0.5);
                Vector3r p1 = p0 - Vector3r(0.0, 0.0, L);
                LinearSpline spline(p0, p1);

                auto beamEnds = engine.createBeam(spline, section, springs, bodyState, beamProps, segments, physics::BeamColliderMode::InterSegmentHull);
                
                engine.addRigidConstraint(beamEnds.first, walkerBody);
                engine.disableCollisionBetween(beamEnds.first, walkerBody);
            }
        }
    }

    void updateScene(physics::PhysicsEngine& /*engine*/, real_t /*t*/) override {
    }
};