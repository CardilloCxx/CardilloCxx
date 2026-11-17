#pragma once

#include "RingGeometry.hpp"
#include "RockfallTest.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace cardillo;

// Drop test on a horizontal rockfall ring net (Geobrugg ROCCO rings). This is the field test of
// Grassl [1] (setup in [3], Fig. 12), which Eberhardt et al. [2] use as benchmark (Sec. 4.2,
// Fig. 19c): a square ring net of 3.9 m x 3.9 m in which every ring is linked with four neighbours is
// shackled to the top of a rigid steel frame, and a block of about 850 kg falls 3 m onto it (gravity -z).
// The shackles are modelled as elastic strings (LinearDistanceConstraint) to a static frame, optionally
// with a rigid shackle around the ring (RockfallTest.hpp).
//
// References
// [1] H. G. Grassl, Experimentelle und numerische Modellierung des dynamischen Trag- und
//     Verformungsverhaltens von hochflexiblen Schutzsystemen gegen Steinschlag, PhD thesis, ETH Zürich, 2002.
// [2] Eberhardt et al., 1D elastoplastic models with isotropic hardening as differential inclusions in
//     primal form, JTCAM, 2026, Sec. 4.2 (supplementary/ring_nets/Eberhardt2026.pdf).
// [3] Z.-H. Zhu, J.-H. Yin, C.-J. Ouyang, D.-Y. Tan, J.-Q. Qin, Modeling a flexible ring net with the
//     discrete element method, J. Eng. Mech. 146(2), doi:10.1061/(ASCE)EM.1943-7889.0001707.
// [4] EOTA, ETAG 027: Guideline for European technical approval of falling rock protection kits,
//     April 2013, Annex A.2, Figure A.2.
// [5] Geobrugg, product data sheet RXE-8000 (BAFU): ROCCO ring net 19/3/300, wire diameter 3 mm, wire
//     tensile strength 1770 N/mm². ROCCO n/d/D denotes rings of n windings of a d mm wire, diameter D mm.
//
// Scene parameters (scene.<name>, defaults in parentheses):
//   rings_x, rings_y (10, 10)       rows of rings_x rings alternate with rows of rings_x - 1 rings;
//                                   rings_y long rows, i.e. 10 x 10 + 9 x 9 = 181 rings
//   net_size_x, net_size_y (3.9)    outer net dimensions in m [1, 3]
//   ring_windings (7), wire_diameter (0.003), ring_inner_diameter (0.3)   ROCCO 7/3/300 [2, 3, 5]
//   ring_segments (16)              beam elements per ring
//   axial_stiffness_factor (0.2)    measured / theoretical axial ring stiffness [1, 3]
//   crossing_clearance (1.5)        gap of linked rings at their crossings in bundle diameters
//   string_length (0.1), string_stiffness (1e7), string_damping (0)   shackle strings (m, N/m, Ns/m)
//   shackle_mode (0)                0: string from the ring node, 1: rigid shackle (RockfallTest.hpp)
//   block_mass (850), block_reference_density (2700), block_mesh (0: rock.obj, 1: rock2.obj)
//   drop_height (3), start_gap (= drop_height)   fall height and initial gap above the net in m
//   extra_boulders (0), extra_size_min (0.15), extra_size_max (0.4), extra_density (2700),
//   extra_spread (2), extra_seed (1)    additional smaller rocks falling with the block
//   primitive_colliders (0)         1: one capsule per ring segment, 0: inter-segment hulls

class NetScene : public SceneBase {
   public:
    const char* sceneName() const override { return "net"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();

        // ---- Rings: ROCCO 7/3/300, 7 windings of a 3 mm wire [2], 300 mm inner diameter, held by three clips [3]
        const real_t windings = cfg.sceneParam("ring_windings", 7);
        const real_t wireRadius = cfg.sceneParam("wire_diameter", 0.003) / 2;
        const real_t innerRadius = cfg.sceneParam("ring_inner_diameter", 0.3) / 2;
        const size_t segments = (size_t)cfg.sceneParam("ring_segments", 16);

        // The windings form a loose bundle. Contact uses its envelope (a hexagonally packed 7-wire bundle
        // is 3 wire diameters wide, fill factor 7/9); mass and stiffness follow the wires: area n A_w and,
        // for bending and torsion, n I_w of separate wires instead of the solid envelope [3, Eq. 10].
        const real_t bundleRadius = 3 * wireRadius * std::sqrt(windings / 7);
        const real_t ringRadius = innerRadius + bundleRadius;  // centreline
        const real_t areaRatio = windings * wireRadius * wireRadius / (bundleRadius * bundleRadius);
        const real_t inertiaRatio = windings * std::pow(wireRadius / bundleRadius, 4);
        // Grassl measured an axial ring stiffness of 1/5 of n E A_w for rings with 7, 12 and 19 windings,
        // as the windings are not loaded evenly [1], [3, Eq. 14].
        const real_t axialFactor = cfg.sceneParam("axial_stiffness_factor", 0.2);
        const real_t E = 2.0e11;  // high-tensile steel wire [3]; tensile strength 1770 N/mm² [5] (no plasticity here)
        const real_t nu = 0.3;
        const real_t steelDensity = 7850;
        const auto section = physics::BeamCrossSection::round(bundleRadius);
        const auto springs = physics::BeamSpringParams::fromMaterial(E, nu, axialFactor * areaRatio, areaRatio, inertiaRatio, inertiaRatio, inertiaRatio);
        const auto ringProps = physics::RigidProps::withDensity(steelDensity * areaRatio);

        // ---- Net [1, 3]: rows of nx and nx - 1 rings alternate, every ring is linked with its four
        // diagonal neighbours (RingNetLattice).
        const RingNetLattice lattice((int)cfg.sceneParam("rings_x", 10), (int)cfg.sceneParam("rings_y", 10), cfg.sceneParam("net_size_x", 3.9), cfg.sceneParam("net_size_y", 3.9), ringRadius,
                                     bundleRadius);

        // Saddle rings (RingGeometry.hpp) start the woven net flat and interlinked; at the crossings the
        // centrelines of linked rings are crossing_clearance bundle diameters apart.
        const real_t clearance = cfg.sceneParam("crossing_clearance", 1.5);
        const real_t saddle = saddleAmplitude(ringRadius, lattice.linkOffset(), clearance * 2 * bundleRadius);
        const SaddleRingSpline ring(ringRadius, saddle);
        const real_t linkGap = polygonDistance(ring, Vector3r::Zero(), ring, lattice.linkOffset(), segments) - 2 * bundleRadius;
        if (linkGap <= 0) throw std::runtime_error("NetScene: linked rings intersect initially; increase scene.crossing_clearance or scene.ring_segments");

        const auto colliderMode = cfg.sceneParam("primitive_colliders", 0) != 0 ? physics::BeamColliderMode::RigidBodyPrimitive : physics::BeamColliderMode::InterSegmentHull;

        const std::vector<RingNetLattice::Site> sites = lattice.sites();
        std::vector<entt::entity> roots;
        for (const RingNetLattice::Site& site : sites) {
            roots.push_back(engine.createBeam(ring, section, springs, physics::RigidState(site.centre, Quaternion4r::Identity()), ringProps, segments, colliderMode).first);
        }

        // ---- Frame and shackles: each boundary ring hangs on the frame, corner rings on both sides, by an
        // elastic string from its outermost node (shackle_mode 0) or by a rigid shackle around it (shackle_mode 1,
        // RockfallTest.hpp).
        const real_t stringLength = cfg.sceneParam("string_length", 0.1);
        const real_t stringStiffness = cfg.sceneParam("string_stiffness", 1e7);
        const real_t stringDamping = cfg.sceneParam("string_damping", 0);
        const int shackleMode = (int)cfg.sceneParam("shackle_mode", 0);
        std::optional<ShackleModel> shackle;
        if (shackleMode == 1) shackle.emplace(cfg, bundleRadius);
        const NetFrame frame(engine, (lattice.nx - 1) * lattice.ax + ringRadius + stringLength, (lattice.ny - 1) * lattice.ay + ringRadius + stringLength);
        size_t strings = 0;
        for (size_t i = 0; i < sites.size(); ++i) {
            for (const Vector3r& d : sites[i].boundaryDirections) {
                if (shackleMode == 1) {
                    const misc::SplineSample s = ring.sample(ringAlpha(d));
                    shackle->add(engine, frame, sites[i].centre + s.position, s.tangent, d, stringLength, stringStiffness, stringDamping);
                } else {
                    const entt::entity node = outermostNode(engine, roots[i], d);
                    frame.addString(engine, node, d, engine.getPosition(node).head<3>() + stringLength * d, stringStiffness, stringDamping);
                }
                ++strings;
            }
        }

        // ---- Block (RockfallTest.hpp)
        addRockfallBlock(engine, saddle + bundleRadius);

        const real_t ringMass = steelDensity * windings * (real_t)M_PI * wireRadius * wireRadius * ring.totalLength();
        std::cout << "Ring net: " << sites.size() << " rings (mass " << ringMass << " kg each, " << (real_t)sites.size() * ringMass << " kg), " << strings << " shackle strings\n"
                  << "  link distance " << lattice.linkDistance() << " m (taut at " << 2 * (ringRadius - bundleRadius) << " m), saddle amplitude " << saddle << " m, link gap " << linkGap << " m\n";
    }
};
