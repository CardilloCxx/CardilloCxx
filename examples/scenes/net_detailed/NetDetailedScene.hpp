#pragma once

#include "../net/WireRing.hpp"
#include <cardillo/physics/assembly/dynamics_assembler.hpp>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cardillo;

// Pull test on five rockfall net rings, every wire of a ring modelled as a Cosserat rod:
//
//   | O     O |
//   |    O    |  ->
//   | O     O |
//
// Two rings hang on pins of each of two rigid blocks (pins between two plates), and a centre ring links all four diagonally (as
// in a 4-neighbour ring net). The left block is fixed, the right block is pulled in +x (no gravity).
// A ring is a bundle of n closed loops of a steel wire: ROCCO 12/3/300 is a 3 mm wire wound 12 times
// to a ring of 300 mm inner diameter, held together by three clips [1, 2]. The wires interact by
// frictional contact. The pull force on the right pins and the reaction force on the left pins over
// the pin displacement go to <output.folder>/<prefix>_force.csv; cf. the ring tests of Grassl [1] in
// [2], Figs. 10 and 11.
//
// Wires and clips: WireRing.hpp (scene.clip_mode: 0 none, 1 pinned, 2 rigid, 3 pre-tensioned band).
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Z.-H. Zhu, J.-H. Yin, C.-J. Ouyang, D.-Y. Tan, J.-Q. Qin, Modeling a flexible ring net with the
//     discrete element method, J. Eng. Mech. 146(2), doi:10.1061/(ASCE)EM.1943-7889.0001707.
//
// Scene parameters (scene.<name>, defaults in parentheses), besides those of WireRing.hpp
// (ring_wire 0, ring_windings 12, wire_segments 48):
//   link_slack (0.01)         initial slack of the links: centre distance below the taut one (m)
//   crossing_clearance (1.5)  gap of linked rings at their crossings in bundle diameters
//   pin_diameter (0.02), pin_gap (0.002)
//   pull_speed (0.5), pull_ramp (0.02)   speed of the right block (m/s), smooth start over pull_ramp (s)

class NetDetailedScene : public SceneBase {
   public:
    const char* sceneName() const override { return "net_detailed"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();

        // ---- Wires and bundle: ROCCO 12/3/300 [1, 2]
        const WireRingModel model(cfg, 0, 12, 48);
        const real_t bundleRadius = model.bundleRadius;
        const real_t ringRadius = model.ringRadius;

        // ---- Rings: centre ring at the origin, corner rings at (±a, ±a), linked with link_slack below
        // the taut centre distance 2 (R - r_b).
        const real_t linkDistance = 2 * (ringRadius - bundleRadius) - cfg.sceneParam("link_slack", 0.01);
        const real_t a = linkDistance / std::sqrt((real_t)2);
        const real_t clearance = cfg.sceneParam("crossing_clearance", 1.5);
        const real_t saddle = saddleAmplitude(ringRadius, Vector3r(a, a, 0), clearance * 2 * bundleRadius);
        const SaddleRingSpline centreline(ringRadius, saddle);
        const WireRingModel::Wires wires = model.wires(centreline);

        // Initial gaps between the wires of a bundle and between the wires of linked rings.
        const real_t bundleGap = model.bundleGap(wires);
        const real_t linkGap = model.linkGap(wires, Vector3r(a, a, 0));
        if (bundleGap <= 0 || linkGap <= 0) throw std::runtime_error("NetDetailedScene: wires intersect initially; increase scene.wire_gap, scene.crossing_clearance or scene.wire_segments");

        const std::array<Vector3r, 5> centres{Vector3r::Zero(), Vector3r(-a, a, 0), Vector3r(-a, -a, 0), Vector3r(a, a, 0), Vector3r(a, -a, 0)};
        for (const Vector3r& centre : centres) model.addRing(engine, centreline, wires, centre);

        // ---- Pins along z through the outer side of the corner rings, each held between two plates
        // above and below the ring (the blocks). Left: fixed, right: pulled in +x.
        const real_t pinRadius = cfg.sceneParam("pin_diameter", 0.02) / 2;
        const real_t pinGap = cfg.sceneParam("pin_gap", 0.002);
        const real_t pinX = a + ringRadius - bundleRadius - pinRadius - pinGap;
        m_pullSpeed = cfg.sceneParam("pull_speed", 0.5);
        m_pullRamp = cfg.sceneParam("pull_ramp", 0.02);
        // The plates cover the bundle around the pin; their faces lie pin_gap above and below the
        // bundle within that footprint (corner ring at (-a, a), pin at (-pinX, a); all corners alike).
        const real_t plateHalfWidth = 2 * (pinRadius + pinGap + bundleRadius);
        real_t zLow = std::numeric_limits<real_t>::infinity();
        real_t zHigh = -std::numeric_limits<real_t>::infinity();
        constexpr int samples = 4096;
        for (int i = 0; i < samples; ++i) {
            const Vector3r q = Vector3r(-a, a, 0) + centreline.sample((real_t)i / samples).position;
            if (std::abs(q.x() + pinX) > plateHalfWidth + bundleRadius || std::abs(q.y() - a) > plateHalfWidth + bundleRadius) continue;
            zLow = std::min(zLow, q.z() - bundleRadius - pinGap);
            zHigh = std::max(zHigh, q.z() + bundleRadius + pinGap);
        }
        const real_t pinZ = (zLow + zHigh) / 2;
        const real_t pinHalfLength = (zHigh - zLow) / 2;
        const Vector3r plateHalfExtents(plateHalfWidth, plateHalfWidth, pinRadius);
        for (const real_t side : {(real_t)-1, (real_t)1}) {
            std::vector<std::pair<entt::entity, Vector3r>> bodies;
            for (const real_t y : {a, -a}) {
                const Vector3r position(side * pinX, y, pinZ);
                const entt::entity pin = engine.addRigidBody(physics::CapsuleShape(pinRadius, pinHalfLength), physics::RigidState(position), physics::RigidProps(0));
                (side > 0 ? m_rightPins : m_leftPins).push_back(pin);
                bodies.emplace_back(pin, position);
                for (const real_t z : {zLow - plateHalfExtents.z(), zHigh + plateHalfExtents.z()}) {
                    const Vector3r platePosition(side * pinX, y, z);
                    bodies.emplace_back(engine.addRigidBody(physics::CubeShape(plateHalfExtents), physics::RigidState(platePosition), physics::RigidProps(0)), platePosition);
                }
            }
            if (side < 0) continue;
            for (const auto& [body, position] : bodies) {
                const real_t v = m_pullSpeed;
                const real_t ramp = m_pullRamp;
                engine.addTrajectory(body, [position, v, ramp](real_t t) -> TrajectoryPose { return TrajectoryPose{position + Vector3r::UnitX() * displacement(v, ramp, t), Quaternion4r::Identity()}; });
            }
        }

        std::filesystem::create_directories(cfg.output_folder);
        m_forceFile.open(cfg.output_folder + "/" + cfg.outputPrefix() + "_force.csv");
        m_forceFile << "t,displacement,pull_force,left_force\n";

        std::cout << "Rings: 5 x " << model.windings << " wires (" << model.ringMass(wires) << " kg per ring), bundle radius " << bundleRadius << " m, gaps: bundle " << bundleGap << " m, link "
                  << linkGap << " m\n"
                  << "  link distance " << linkDistance << " m, saddle amplitude " << saddle << " m, clip mode " << model.clipMode << "\n";
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t) override {
        // Contact forces of the last step on the pins (the impulse acts on body b of a contact).
        Vector3r right = Vector3r::Zero();
        Vector3r left = Vector3r::Zero();
        for (const auto& c : engine.dynamicsAssembler().contacts()) {
            const Vector3r impulse = c.last_impulse.x() * c.normal + c.last_impulse.y() * c.tangent1 + c.last_impulse.z() * c.tangent2;
            auto add = [&](const std::vector<entt::entity>& pins, Vector3r& sum) {
                if (std::find(pins.begin(), pins.end(), c.b) != pins.end()) sum += impulse;
                if (std::find(pins.begin(), pins.end(), c.a) != pins.end()) sum -= impulse;
            };
            add(m_rightPins, right);
            add(m_leftPins, left);
        }
        const real_t h = engine.timeStep();
        m_forceFile << t << ',' << displacement(m_pullSpeed, m_pullRamp, t) << ',' << -right.x() / h << ',' << left.x() / h << '\n';
    }

   private:
    // Displacement of the right block: speed v reached smoothly (constant acceleration) over the ramp time.
    static real_t displacement(real_t v, real_t ramp, real_t t) { return t < ramp ? v * t * t / (2 * ramp) : v * (t - ramp / 2); }

    std::vector<entt::entity> m_rightPins;
    std::vector<entt::entity> m_leftPins;
    real_t m_pullSpeed{(real_t)0};
    real_t m_pullRamp{(real_t)0};
    std::ofstream m_forceFile;
};
