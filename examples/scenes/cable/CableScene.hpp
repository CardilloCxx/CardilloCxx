#pragma once

#include <Eigen/Geometry>
#include <cmath>
#include <algorithm>
#include <optional>
#include <vector>
#include "../SceneBase.hpp"
#include "misc/spline.hpp"

using namespace cardillo;

class CableScene : public SceneBase {
public:
    const char* sceneName() const override { return "cable"; }
    CableScene() = default;
    ~CableScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Set standard gravity (Z is up)
        engine.setGravity(Vector3r(0,  -9.81, 0));

        // Add a ground plane somewhat below the jack to catch the longer cable
        physics::CubeShape groundShape(Vector3r(10.0, 0.1, 10.0));
        engine.addStaticBody(groundShape, physics::RigidState(Vector3r(0.0, -0.5, 0.0)));

        physics::MeshShape jackShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/ethernet_jack.obj");
        engine.addStaticBody(jackShape, physics::RigidState(Vector3r::Zero()));

        Vector3r loosePortPos(-0.5, 0.0, 0.3);
        Vector3r initialPlugPos(0.0, 0.0, 0.6);

        std::vector<Vector3r> cablePoints = {
            loosePortPos,
            Vector3r(0.0, 0.2, 0.4),
            Vector3r(0.3, -0.1, 0.3),
            Vector3r(-0.1, 0.2, 0.1),
            initialPlugPos
        };
        
        CatmullRomSpline cableSpline(cablePoints, false);

        // --------------------------------------------------------
        physics::MeshShape plugShape(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/ethernet_plug.obj");
        physics::RigidProps movingPlugProps;
        movingPlugProps.mass = (real_t)0.0; 
        movingPlugProps.collidable = false; 
        
        SplineSample cableEndSample = cableSpline.sample(1.0);
        
        // The cable enters from -Z in untransformed space. 
        Quaternion4r initialPlugRot = Eigen::Quaternion<real_t>::FromTwoVectors(Vector3r::UnitZ(), -cableEndSample.tangent);
        entt::entity movingPlug = engine.addRigidBody(plugShape, physics::RigidState(initialPlugPos, initialPlugRot), movingPlugProps);

        // Define the insertion trajectory via a spline
        std::vector<Vector3r> trajPoints = {
            initialPlugPos,
            Vector3r(0.15, -0.1, 0.4), // Swoop to the side
            Vector3r(0.0, 0.0, 0.15),  // Align straight above the jack
            Vector3r::Zero()           // Final docked position
        };
        CatmullRomSpline trajSpline(trajPoints, false);

        std::optional<std::function<TrajectoryPose(real_t)>> poseFunc = [trajSpline, initialPlugRot](real_t t) -> TrajectoryPose {
            real_t time_clamped = std::min(t, (real_t)4.0);
            real_t progress = time_clamped / (real_t)4.0;
            
            real_t s = progress * progress * (3 - 2 * progress);
            
            TrajectoryPose pose;
            pose.first = trajSpline.sample(s).position; 
            pose.second = initialPlugRot.slerp(s, Quaternion4r::Identity());
            
            return pose;
        };
        engine.addTrajectory(movingPlug, poseFunc, std::nullopt);

        physics::RigidProps loosePortProps = physics::RigidProps::withDensity((real_t)500.0);
        loosePortProps.collidable = false; 
        
        // Sample the beginning of the cable spline (alpha = 0.0
        SplineSample cableStartSample = cableSpline.sample(0.0);
        
        // At the start of the spline, the cable runs in the +tangent direction.
        Quaternion4r looseRot = Eigen::Quaternion<real_t>::FromTwoVectors(Vector3r::UnitZ(), cableStartSample.tangent);
        entt::entity loosePort = engine.addRigidBody(plugShape, physics::RigidState(loosePortPos, looseRot), loosePortProps);

        const real_t r_cable = (real_t)0.003;
        const real_t d = (real_t)(2 * r_cable);
        
        physics::BeamCrossSection cableSection(d, d, physics::BeamBodyType::Capsule);
        auto cableSprings = physics::BeamSpringParams::fromMaterial((real_t)5e7, (real_t)0.4);
        cableSprings.setDampingFromFactor((real_t)1.0);
        cableSprings.kappa0 = Vector3r::Zero();   // Set relaxed curvature to zero for a straight cable
        
        physics::RigidProps cableProps = physics::RigidProps::withDensity((real_t)1200.0);
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());

        const size_t cableSegments = 320;
        auto cableEnds = engine.createBeam(cableSpline, cableSection, cableSprings, stateDefaults, cableProps, cableSegments);
        
        // First end attached to the loose dynamic port
        engine.addRigidConstraint(cableEnds.first, loosePort);
        // Second end attached to the moving, trajectory-driven plug
        engine.addRigidConstraint(cableEnds.second, movingPlug);
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t t, real_t dt) override {
        (void)engine;
        (void)t;
        (void)dt;
    }
};