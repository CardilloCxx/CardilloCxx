#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace cardillo;

// Ciliary walkers on a vertically vibrated plate after S. Mohanty, P. Baconnier, H. A. H. Schomaker,
// A. Comoretto, M. van Hecke, J. T. B. Overvelde, "Multimodal motion and behavior switching of
// multistable ciliary walkers", arXiv:2512.09840 (2025). Dimensions from Fig. S1, material from
// Table SI; a silicone (ED22F) walker of ~9 g stands on flexible cilia that buckle under its weight
// (and an added top weight), and the buckled state sets the direction of propulsion.
//
// Variants (`scene.walker`):
//   bidirectional  28 x 28 x 10 mm body, 5 x 5 rectangular cilia (6 mm, width 1.2 -> 0.2 mm,
//                  depth 1.5 mm, pitch 6 mm), soft (buckling) direction along e_par = x (Fig. S1A)
//   isotropic      same body, 5 x 5 cylindrical cilia (5 mm: diameter 1.0 mm over 3 mm, then
//                  0.6 -> 0.5 mm), pitch 6 mm (Fig. S1B); grid_n = 10, pitch = 3 mm gives the
//                  10 x 10 walker of Fig. 5C
//   ring           disk-shaped body with 16 rectangular cilia on a ring, soft direction along e_theta
//                  (Fig. 5B, S11B); body diameter and ring radius are not given in the paper and are
//                  assumptions (scene parameters)
// The plate is the 16 cm aluminium disk of Fig. S2, the arena wall a ring of the same diameter that
// moves with the plate. A small tip offset along the soft direction selects the buckling direction.
// `scene.n_walkers` > 1 places several identical walkers in the arena (swarm experiments, Fig. 6).
class CiliaWalkerScene : public SceneBase {
   public:
    const char* sceneName() const override { return "cilia_walker"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();
        const auto param = [&](const char* name, real_t fallback) { return cfg.sceneParam(name, fallback); };
        const auto walkerIt = cfg.scene_params.find("walker");
        const std::string walker = walkerIt == cfg.scene_params.end() ? "bidirectional" : walkerIt->second;
        if (walker != "bidirectional" && walker != "isotropic" && walker != "ring")
            throw std::runtime_error("CiliaWalkerScene: scene.walker must be bidirectional, isotropic or ring, got '" + walker + "'");
        const bool isotropic = walker == "isotropic";
        const bool ring = walker == "ring";

        engine.setGravity(Vector3r(0, 0, -9.81));

        // --- Material (ED22F silicone; E from Table SI, density from the 9 g walker of Fig. S1A)
        const real_t E = param("youngs_modulus", (real_t)1.2e6);
        const real_t nu = param("poisson_ratio", (real_t)0.49);
        const real_t rho = param("density", (real_t)1120);

        auto springs = physics::BeamSpringParams::fromMaterial(E, nu);
        // Stiffness-proportional material damping of the rod elements: D = d * K for every element
        // row, d in seconds (modal damping ratio zeta = d * omega / 2).
        springs.setDampingFromFactor(param("rod_damping_factor", (real_t)0));
        const physics::RigidProps ciliumProps = physics::RigidProps::withDensity(rho);

        // --- Cilium geometry
        const real_t L = param("cilium_length", isotropic ? (real_t)5e-3 : (real_t)6e-3);
        const size_t segments = (size_t)param("cilium_segments", isotropic ? (real_t)10 : (real_t)12);
        const auto section = isotropic ? physics::BeamCrossSection::procedural([](real_t s) {
            // diameter 1.0 mm over the first 3 mm, then tapering from 0.6 to 0.5 mm (L = 5 mm)
            const real_t z = s * (real_t)5e-3;
            const real_t d = z < (real_t)3e-3 ? (real_t)1.0e-3 : (real_t)0.6e-3 - (real_t)0.1e-3 * (z - (real_t)3e-3) / (real_t)2e-3;
            return physics::BeamCrossSection::round((real_t)0.5 * d, 12);
        })
        : physics::BeamCrossSection::procedural([](real_t s) {
            // width (soft bending direction) 1.2 -> 0.2 mm, depth 1.5 mm
            return physics::BeamCrossSection::rectangular((real_t)1.2e-3 - (real_t)1.0e-3 * s, (real_t)1.5e-3);
        });
        // Lateral tip offset that selects the buckling direction: along the soft direction (e_par
        // for the bidirectional, e_theta for the ring walker), along tip_offset_angle for isotropic cilia.
        const real_t tipOffset = param("tip_offset", (real_t)1e-4);
        const real_t tipOffsetAngle = param("tip_offset_angle_deg", (real_t)0) * (real_t)M_PI / (real_t)180;

        // --- Vibrating plate and arena (Fig. S2: 16 cm aluminium disk)
        const real_t frequency = param("frequency", (real_t)40);
        const real_t gamma = param("gamma", (real_t)2.45);
        const real_t amplitude = gamma * (real_t)9.81 / std::pow((real_t)2 * (real_t)M_PI * frequency, 2);
        const auto lift = [=](real_t t) { return amplitude * std::sin((real_t)2 * (real_t)M_PI * frequency * t); };

        const real_t plateRadius = (real_t)0.5 * param("plate_diameter", (real_t)0.16);
        const real_t plateThickness = (real_t)0.01;
        physics::RigidProps drivenProps;
        drivenProps.mass = std::nullopt;
        drivenProps.visual = true;
        const Vector3r plateCenter(0, 0, -(real_t)0.5 * plateThickness);  // plate top at z = 0
        const entt::entity plate = engine.addRigidBody(physics::CylinderShape(plateRadius, (real_t)0.5 * plateThickness), physics::RigidState(plateCenter), drivenProps);
        engine.addTrajectory(plate, [=](real_t t) { return TrajectoryPose{plateCenter + Vector3r(0, 0, lift(t)), Quaternion4r::Identity()}; });

        const real_t arenaRadius = (real_t)0.5 * param("arena_diameter", (real_t)0.16);
        const real_t wallHeight = (real_t)0.03;
        const real_t wallThickness = (real_t)0.008;
        const int numWallSegments = 32;
        const real_t wallSegmentLength = (real_t)2 * (arenaRadius + wallThickness) * std::tan((real_t)M_PI / numWallSegments);
        for (int i = 0; i < numWallSegments; ++i) {
            const real_t phi = ((real_t)i + (real_t)0.5) / numWallSegments * (real_t)2 * (real_t)M_PI;
            const Vector3r wallCenter((arenaRadius + (real_t)0.5 * wallThickness) * std::cos(phi), (arenaRadius + (real_t)0.5 * wallThickness) * std::sin(phi),
                                      (real_t)0.5 * wallHeight);
            const Quaternion4r wallRot(Eigen::AngleAxis<real_t>(phi, Vector3r::UnitZ()));
            physics::CubeShape wallShape(Vector3r((real_t)0.5 * wallThickness, (real_t)0.5 * wallSegmentLength, (real_t)0.5 * wallHeight));
            const entt::entity wall = engine.addRigidBody(wallShape, physics::RigidState(wallCenter, Vector3r::Zero(), wallRot), drivenProps);
            engine.addTrajectory(wall, [=](real_t t) { return TrajectoryPose{wallCenter + Vector3r(0, 0, lift(t)), wallRot}; });
            engine.disableCollisionBetween(wall, plate);
        }

        // --- Walkers: body (silicone), added top weight and cilia. Several walkers are placed on a
        // circle of radius walker_placement_radius around the arena center, walker k rotated by
        // k * walker_yaw_step_deg about the vertical axis.
        const int numWalkers = std::max(1, (int)param("n_walkers", (real_t)1));
        const real_t placementRadius = numWalkers > 1 ? param("walker_placement_radius", (real_t)0.04) : (real_t)0;
        const real_t yawStep = param("walker_yaw_step_deg", (real_t)90) * (real_t)M_PI / (real_t)180;

        const real_t bodyHeight = (real_t)0.01;
        const real_t bodyWidth = (real_t)0.028;  // square walkers
        const real_t bodyRadius = (real_t)0.5 * param("body_diameter", (real_t)0.035);  // ring walker (assumed)
        const real_t ciliumHeight = std::sqrt(L * L - tipOffset * tipOffset);
        const real_t bodyZ = ciliumHeight + (real_t)1e-5 + (real_t)0.5 * bodyHeight;  // cilia tips just above the plate
        const physics::RigidProps bodyProps = physics::RigidProps::withDensity(rho);
        const real_t addedWeight = param("added_weight", (real_t)0.012);

        // Cilia: base positions (body frame) and rotation about the cilium axis that turns the soft
        // direction (world x for an unrotated vertical cilium) into the desired direction.
        struct CiliumPlacement {
            Vector2r base;
            real_t rotation;
        };
        std::vector<CiliumPlacement> placements;
        if (ring) {
            const int n = (int)param("n_cilia", (real_t)16);
            const real_t ringRadius = param("ring_radius", (real_t)0.0145);
            for (int i = 0; i < n; ++i) {
                const real_t phi = (real_t)i / n * (real_t)2 * (real_t)M_PI;
                placements.push_back({ringRadius * Vector2r(std::cos(phi), std::sin(phi)), phi + (real_t)0.5 * (real_t)M_PI});  // x -> e_theta
            }
        } else {
            const int n = (int)param("grid_n", (real_t)5);
            const real_t pitch = param("pitch", (real_t)0.006);
            const real_t offset = (real_t)0.5 * (n - 1) * pitch;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) placements.push_back({Vector2r(-offset + i * pitch, -offset + j * pitch), (real_t)0});
        }
        const Vector3r tipShift = tipOffset * Vector3r(std::cos(isotropic ? tipOffsetAngle : (real_t)0), std::sin(isotropic ? tipOffsetAngle : (real_t)0), 0);

        for (int k = 0; k < numWalkers; ++k) {
            const real_t angle = (real_t)2 * (real_t)M_PI * k / numWalkers;
            const real_t yaw = k * yawStep;
            const Matrix33r Ryaw = Eigen::AngleAxis<real_t>(yaw, Vector3r::UnitZ()).toRotationMatrix();
            const Quaternion4r qyaw(Ryaw);
            const Vector3r bodyPos(placementRadius * std::cos(angle), placementRadius * std::sin(angle), bodyZ);
            const physics::RigidState bodyState(bodyPos, Vector3r::Zero(), qyaw);
            const entt::entity body = ring ? engine.addRigidBody(physics::CylinderShape(bodyRadius, (real_t)0.5 * bodyHeight), bodyState, bodyProps)
                                           : engine.addRigidBody(physics::CubeShape(Vector3r((real_t)0.5 * bodyWidth, (real_t)0.5 * bodyWidth, (real_t)0.5 * bodyHeight)), bodyState, bodyProps);
            // position/velocity/orientation in <output.folder>/cilia_walker_tracked.csv
            engine.track(body, numWalkers == 1 ? std::string("walker") : "walker_" + std::to_string(k));

            // Added top weight (VeroClear, press-fit in the experiments): same footprint as the body
            if (addedWeight > (real_t)0) {
                const real_t footprint = ring ? (real_t)M_PI * bodyRadius * bodyRadius : bodyWidth * bodyWidth;
                const real_t weightHeight = addedWeight / ((real_t)1180 * footprint);
                physics::RigidProps weightProps;
                weightProps.mass = addedWeight;
                const physics::RigidState weightState(bodyPos + Vector3r(0, 0, (real_t)0.5 * (bodyHeight + weightHeight)), Vector3r::Zero(), qyaw);
                const entt::entity weight = ring ? engine.addRigidBody(physics::CylinderShape(bodyRadius, (real_t)0.5 * weightHeight), weightState, weightProps)
                                                 : engine.addRigidBody(physics::CubeShape(Vector3r((real_t)0.5 * bodyWidth, (real_t)0.5 * bodyWidth, (real_t)0.5 * weightHeight)), weightState, weightProps);
                engine.addRigidConstraint(weight, body);
                engine.disableCollisionBetween(weight, body);
            }

            for (const auto& p : placements) {
                const Vector2r base = Ryaw.topLeftCorner<2, 2>() * p.base;
                const Vector3r p0(base.x(), base.y(), -(real_t)0.5 * bodyHeight);
                const Vector3r p1 = p0 - Vector3r(0, 0, ciliumHeight) + tipShift;
                // createBeam rotates the cilium about its own center by the orientation of this state,
                // which also turns the tip offset into the soft direction of the yawed walker
                const physics::RigidState ciliumState(bodyPos, Vector3r::Zero(), Quaternion4r(Eigen::AngleAxis<real_t>(yaw + p.rotation, Vector3r::UnitZ())));
                const auto ends = engine.createBeam(misc::LinearSpline(p0, p1), section, springs, ciliumState, ciliumProps, segments, physics::BeamColliderMode::InterSegmentHull);
                engine.addRigidConstraint(ends.first, body);
                // The cilium is attached to the body: no contacts between any of its segments and the
                // body (the hulls next to the root touch the bottom face of the body).
                for (std::optional<entt::entity> seg = ends.first; seg; seg = engine.ecs().get<C_BeamElement>(*seg).next) engine.disableCollisionBetween(*seg, body);
            }
            if (k == 0)
                std::cout << "Ciliary walker '" << walker << "': " << numWalkers << " x (" << placements.size() << " cilia, body " << engine.getMass(body)(0, 0) * 1e3
                          << " g + added weight " << addedWeight * 1e3 << " g), plate amplitude " << amplitude * 1e3 << " mm" << std::endl;
        }
    }
};
