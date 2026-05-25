#pragma once

#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include "../SceneBase.hpp"
#include "misc/spline.hpp"

using namespace cardillo;

class CableBunchScene : public SceneBase {
   public:
    const char* sceneName() const override { return "cable_bunch"; }
    CableBunchScene() = default;
    ~CableBunchScene() = default;

    void populate(cardillo::physics::PhysicsEngine& engine) override {
        using namespace cardillo;
        using namespace cardillo::misc;

        // Gravity on (Z is up)
        engine.setGravity(Vector3r(0, 0, -9.81));

        // Create a central arbitrary spline for the cable bunch shape
        std::vector<Vector3r> centerPoints = {Vector3r(0.0, 0.0, 0.05), Vector3r(0.2, 0.0, 0.05), Vector3r(0.4, 0.0, 0.05),
                                              Vector3r(0.6, 0.0, 0.05), Vector3r(0.8, 0.0, 0.05), Vector3r(1.0, 0.0, 0.05)};
        CatmullRomSpline centerSpline(centerPoints, false);

        // Beam / cable physical properties
        const real_t r_cable = (real_t)0.005;  // 5mm radius
        const real_t d = (real_t)(2 * r_cable);
        const real_t E = (real_t)1e7;  // Softer material for cables
        const real_t nu = (real_t)0.30;
        const real_t rho = (real_t)1500.0;

        physics::BeamCrossSection section(d, d, physics::BeamBodyType::Capsule);
        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        springs.setDampingFromFactor((real_t)1.0);
        physics::RigidProps props = physics::RigidProps::withDensity(rho);
        physics::RigidState stateDefaults(Vector3r::Zero(), Vector3r::Zero(), Quaternion4r::Identity());

        const size_t segments = 50;

        // Bunch setup
        const int numCables = 6;
        const real_t bunchRadius = r_cable * (real_t)2.5;

        // Generating offset splines
        std::vector<CatmullRomSpline> allSplines;
        allSplines.push_back(centerSpline);  // The core wire

        const int numSamples = 25;
        for (int i = 0; i < numCables; ++i) {
            real_t theta = (real_t)i / (real_t)numCables * (real_t)(2.0 * M_PI);
            std::vector<Vector3r> offsetPoints;

            for (int s = 0; s < numSamples; ++s) {
                real_t alpha = (real_t)s / (real_t)(numSamples - 1);
                SplineSample samp = centerSpline.sample(alpha);

                // Offset in normal and binormal directions
                Vector3r offsetPos = samp.position + samp.normal * (bunchRadius * std::cos(theta)) + samp.binormal * (bunchRadius * std::sin(theta));

                offsetPoints.push_back(offsetPos);
            }
            allSplines.emplace_back(offsetPoints, false);
        }

        // Invisible tiny anchor used for the trajectory (no visual, no collisions)
        physics::CubeShape anchorShape(Vector3r((real_t)0.001, (real_t)0.001, (real_t)0.001));
        Vector3r wallPos = centerPoints.back() + Vector3r(0.02, 0, 0);

        physics::RigidProps anchorProps;
        anchorProps.mass = (real_t)0.0;  // static/kinematic
        anchorProps.visual = false;      // do not render
        anchorProps.collidable = false;  // avoid collisions for the tiny proxy
        entt::entity anchor = engine.addRigidBody(anchorShape, physics::RigidState(wallPos), anchorProps);
        m_anchor = anchor;  // Save it to apply trajectories later
        m_startPos = wallPos;

        // Visible cube representing the bundle that follows the anchor
        const real_t visDiam = (bunchRadius * (real_t)2.0) * (real_t)1.4;
        physics::CubeShape visShape(Vector3r(visDiam, visDiam, visDiam));
        physics::RigidProps visProps = physics::RigidProps::withDensity((real_t)100.0);
        entt::entity visCube = engine.addRigidBody(visShape, physics::RigidState(centerPoints.front()), visProps);

        // Create physical beams in the engine and attach their ends to the (invisible) anchor
        for (const auto& spl : allSplines) {
            auto ends = engine.createBeam(spl, section, springs, stateDefaults, props, segments);
            engine.addRigidConstraint(ends.second, anchor);
            engine.addRigidConstraint(ends.first, visCube);
        }

        // Add cable ties (rings with flat rectangular cross section) along the center spline
        const real_t tieSpacing = 0.2;  // Tie every 20cm
        const int numTies = std::floor(centerSpline.totalLength() / tieSpacing);
        // Belt cross section for the zip-ties
        const real_t w_tie = 0.001;  // Radial thickness
        const real_t h_tie = 0.010;  // Logitudinal width
        physics::BeamCrossSection tieSection(w_tie, h_tie, physics::BeamBodyType::Cube);
        physics::RigidProps tieProps = physics::RigidProps::withDensity(1000.0);

        auto tieSprings = physics::BeamSpringParams::fromMaterial(1e9, 0.4);
        tieSprings.gamma0 = Vector3r::Zero();
        tieSprings.gamma0->x() = -h_tie * 0.05;  // Slight pre-compression to keep the ties tight around the bundle
        tieSprings.kappa0 = Vector3r::Zero();

        for (int i = 1; i < numTies; ++i) {
            real_t alpha = (real_t)i / (real_t)numTies;
            SplineSample samp = centerSpline.sample(alpha);

            // Generate a circular spline around the bundle
            // Spline lies in the normal-binormal plane.
            CircleSpline circleSpline(samp.position, bunchRadius * 1.4, samp.tangent, samp.normal);

            // Loop = true for creating a continuous ring-belt
            engine.createBeam(circleSpline, tieSection, tieSprings, stateDefaults, tieProps, 20);
        }

        // Obstacles (cylinders)
        std::vector<Vector3r> obsPos = {Vector3r(1.5, 0.1, 0.05), Vector3r(2.0, -0.1, 0.05), Vector3r(2.5, 0.1, 0.05)};
        const real_t obsR = 0.05;
        const real_t obsH = 0.2;
        physics::CylinderShape cylShape(obsR, obsH);
        for (auto& p : obsPos) {
            engine.addStaticBody(cylShape, physics::RigidState(p));
        }

        // Add Floor
        engine.addStaticBody(physics::CubeShape(Vector3r(10.0, 10.0, 0.05)), physics::RigidState(Vector3r(0.0, 0.0, -0.05)));

        // Trajectory spline through the obstacles (visualized by a static, non-collidable beam)
        std::vector<Vector3r> trajPoints = {wallPos, Vector3r(1.5, 0.0, 0.05), Vector3r(2.0, -0.15, 0.05), Vector3r(2.5, 0.15, 0.05), Vector3r(3.0, 0.0, 0.05)};
        CatmullRomSpline trajSpline(trajPoints, false);

        physics::BeamCrossSection visSection((real_t)0.01, (real_t)0.01, physics::BeamBodyType::Cube);
        auto visSprings = physics::BeamSpringParams::fromMaterial((real_t)1e8, (real_t)0.3);
        physics::RigidProps visBeamProps;  // default: no mass/density -> static
        visBeamProps.collidable = false;   // purely visual
        visBeamProps.visual = true;
        engine.createBeam(trajSpline, visSection, visSprings, stateDefaults, visBeamProps, 100);

        // Drive the invisible anchor along the trajectory spline with position only.
        engine.addTrajectory(anchor, trajSpline, (real_t)6.0);
    }

    void updateScene(cardillo::physics::PhysicsEngine& engine, real_t /*t*/, real_t /*dt*/) override { (void)engine; }

   private:
    entt::entity m_anchor;
    Vector3r m_startPos;
};
