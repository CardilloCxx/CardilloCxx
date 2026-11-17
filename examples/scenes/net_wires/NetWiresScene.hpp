#pragma once

#include "../net/RockfallTest.hpp"
#include "../net/WireRing.hpp"
#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cardillo;

// Drop test on a horizontal rockfall ring net whose rings are modelled wire by wire: the setup of
// NetScene (Grassl's field test [1], benchmark of Eberhardt et al. [2]: a ring net in which every ring is
// linked with four neighbours is shackled to a rigid frame, a block of about 850 kg falls 3 m onto it,
// gravity -z), with every ring a bundle of wires held together by clips (WireRing.hpp). Default: one
// 3 mm wire wound 9 times to 300 mm inner diameter, three pinned clip sleeves, on the 3.9 m x 3.9 m net
// of [1] with 10 x 10 + 9 x 9 rings.
//
// Net, wires and clips, shackles and block are set by scene parameters, so the scene scales from a few
// rings to large nets (the default net with 32 elements per turn has 52 309 wire nodes).
//
// Scene parameters (scene.<name>, defaults in parentheses), besides those of WireRing.hpp
// (ring_wire 1, ring_windings 9, wire_segments 32, clip_mode, clip_count, band_*, ...) and RockfallTest.hpp
// (block_*, drop_height, start_gap, extra_*, shackle_*):
//   rings_x, rings_y (10, 10)       rows of rings_x rings alternate with rows of rings_x - 1 rings;
//                                   rings_y long rows
//   net_size_x, net_size_y (3.9)    outer net dimensions (m); the frame (mounting) lies string_length
//                                   outside the outermost ring centrelines
//   crossing_clearance (1.5)        gap of linked rings at their crossings in bundle diameters
//   string_length (0.1), string_stiffness (1e7), string_damping (0)   shackle strings (m, N/m, Ns/m)
//   shackle_mode (1)                0: a string from the outermost node of every wire,
//                                   1: rigid shackle around the bundle (RockfallTest.hpp)
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Eberhardt et al., 1D elastoplastic models with isotropic hardening as differential inclusions in
//     primal form, JTCAM, 2026, Sec. 4.2 (supplementary/ring_nets/Eberhardt2026.pdf).

class NetWiresScene : public SceneBase {
   public:
    const char* sceneName() const override { return "net_wires"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();

        // ---- Rings and net
        const WireRingModel model(cfg, 1, 9, 32);
        const RingNetLattice lattice((int)cfg.sceneParam("rings_x", 10), (int)cfg.sceneParam("rings_y", 10), cfg.sceneParam("net_size_x", 3.9), cfg.sceneParam("net_size_y", 3.9), model.ringRadius,
                                     model.bundleRadius);
        const real_t clearance = cfg.sceneParam("crossing_clearance", 1.5);
        const real_t saddle = saddleAmplitude(model.ringRadius, lattice.linkOffset(), clearance * 2 * model.bundleRadius);
        const SaddleRingSpline centreline(model.ringRadius, saddle);
        const WireRingModel::Wires wires = model.wires(centreline);

        // Initial gaps: wires of a bundle, wires of linked rings (all links are alike for identical saddle
        // rings), clips and the wires of the linked rings.
        const std::vector<Vector3r> neighbours{Vector3r(lattice.ax, lattice.ay, 0), Vector3r(-lattice.ax, lattice.ay, 0), Vector3r(lattice.ax, -lattice.ay, 0), Vector3r(-lattice.ax, -lattice.ay, 0)};
        const real_t bundleGap = model.bundleGap(wires);
        const real_t linkGap = model.linkGap(wires, lattice.linkOffset());
        const real_t clipGap = model.clipGap(centreline, wires, neighbours);
        if (bundleGap <= 0 || linkGap <= 0) throw std::runtime_error("NetWiresScene: wires intersect initially; increase scene.wire_gap, scene.crossing_clearance or scene.wire_segments");
        if (clipGap <= 0) throw std::runtime_error("NetWiresScene: clips intersect linked rings initially; change scene.clip_offset_deg");

        const std::vector<RingNetLattice::Site> sites = lattice.sites();
        std::vector<WireRingModel::RingNodes> rings;
        for (const RingNetLattice::Site& site : sites) rings.push_back(model.addRing(engine, centreline, wires, site.centre));

        // ---- Frame and shackles: each boundary ring hangs on the frame, corner rings on both sides.
        const real_t stringLength = cfg.sceneParam("string_length", 0.1);
        const real_t stringStiffness = cfg.sceneParam("string_stiffness", 1e7);
        const real_t stringDamping = cfg.sceneParam("string_damping", 0);
        const int shackleMode = (int)cfg.sceneParam("shackle_mode", 1);
        std::optional<ShackleModel> shackle;
        if (shackleMode == 1) shackle.emplace(cfg, model.bundleRadius);
        if (shackle && model.clipMode != 0) {
            // A shackle and a clip around the same bundle must not overlap along it.
            for (const Vector3r& d : std::array<Vector3r, 4>{Vector3r::UnitX(), Vector3r::UnitY(), -Vector3r::UnitX(), -Vector3r::UnitY()}) {
                for (int c = 0; c < model.clipCount; ++c) {
                    const real_t distance = (centreline.sample(ringAlpha(d)).position - centreline.sample(model.nodeAlpha(model.clipNode(c))).position).norm();
                    if (distance <= model.bandWidth / 2 + shackle->rodRadius()) throw std::runtime_error("NetWiresScene: a clip overlaps a shackle; change scene.clip_offset_deg");
                }
            }
        }
        const NetFrame frame(engine, (lattice.nx - 1) * lattice.ax + model.ringRadius + stringLength, (lattice.ny - 1) * lattice.ay + model.ringRadius + stringLength);
        size_t shackles = 0;
        size_t strings = 0;
        for (size_t i = 0; i < sites.size(); ++i) {
            for (const Vector3r& d : sites[i].boundaryDirections) {
                const misc::SplineSample s = centreline.sample(ringAlpha(d));
                const Vector3r p = sites[i].centre + s.position;
                if (shackle) {
                    shackle->add(engine, frame, p, s.tangent, d, stringLength, stringStiffness, stringDamping);
                    ++strings;
                } else {
                    for (const entt::entity node : model.nodesAt(rings[i], model.nearestNode(ringAlpha(d)))) {
                        frame.addString(engine, node, d, p + stringLength * d, stringStiffness, stringDamping);
                        ++strings;
                    }
                }
                ++shackles;
            }
        }

        // ---- Block (RockfallTest.hpp)
        addRockfallBlock(engine, saddle + model.bundleRadius);

        // ---- Model size: rigid bodies (7 position and 6 velocity coordinates each) and bilateral rows
        // (6 per rod element, 3 per pinned and 6 per rigid clip node, 1 per string); all rings are alike.
        const size_t ringCount = sites.size();
        size_t wireNodes = 0;
        for (const auto& wireNodesOfRing : rings.front()) wireNodes += ringCount * wireNodesOfRing.size();
        const size_t wireElements = ringCount * wires.size() * model.wireSegments();
        size_t clipBodies = 0;
        size_t clipRows = 0;
        if (model.clipMode == 3) clipBodies = ringCount * (size_t)model.clipCount * model.bandSegments;
        if (model.clipMode == 1 || model.clipMode == 2) {
            clipBodies = ringCount * (size_t)model.clipCount;
            for (int c = 0; c < model.clipCount; ++c) clipRows += ringCount * (model.clipMode == 1 ? 3 : 6) * model.nodesAt(rings.front(), model.clipNode(c)).size();
        }
        const size_t shackleBodies = shackle ? shackles : 0;
        const size_t bodies = wireNodes + clipBodies + shackleBodies;
        const size_t bilateralRows = 6 * (wireElements + (model.clipMode == 3 ? clipBodies : 0)) + clipRows + strings;
        std::cout << "Ring net: " << ringCount << " rings, " << model.windings << " windings each " << (model.wound ? "(one wound wire)" : "(closed loops)") << ", " << model.ringMass(wires)
                  << " kg per ring, " << (real_t)ringCount * model.ringMass(wires) << " kg; " << shackles << " shackles, clip mode " << model.clipMode << "\n"
                  << "  bundle radius " << model.bundleRadius << " m, link distance " << lattice.linkDistance() << " m (taut at " << 2 * (model.ringRadius - model.bundleRadius)
                  << " m), saddle amplitude " << saddle << " m\n"
                  << "  initial gaps: bundle " << bundleGap << " m, link " << linkGap << " m" << (model.clipMode != 0 ? ", clip " + std::to_string(clipGap) + " m" : std::string()) << "\n"
                  << "  rigid bodies: " << wireNodes << " wire nodes + " << clipBodies << (model.clipMode == 3 ? " band nodes" : " sleeves") << " + " << shackleBodies << " shackles; positions "
                  << 7 * bodies << ", velocities " << 6 * bodies << ", bilateral rows " << bilateralRows << " (without block and debris)\n";
    }
};
