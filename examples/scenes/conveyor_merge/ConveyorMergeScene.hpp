#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>

using namespace cardillo;

// Infeed-merge conveyor transfer with a waterfall (simplified flat junction geometry): a single
// parcel moves on the infeed belt, which joins the merge belt at an angle of 30 degrees (top view),
// leaves the infeed over its edge along the merge belt, falls by the waterfall height onto the
// merge belt and is carried away. The quantity of interest is the rotation (yaw) of the parcel
// caused by the transfer and the impact.
//
// Coordinates: z up; the merge belt runs along -x with speed v_merge, its surface at z = 0, width
// 1.3 m centered at y = 0. The infeed lies at y > 0.65 m, its surface at z = h (waterfall height),
// and moves along d = (-cos(alpha), -sin(alpha), 0) with v_infeed = v_merge / cos(alpha), so that
// both belts have the same velocity component along the merge direction.
//
// Belts are static bodies with a prescribed surface velocity (engine.setSurfaceVelocity()): the
// geometry is at rest, only the belt surface moves, which enters the contact laws. As in the
// schematic, the infeed is a strip of width 0.8 m at the angle alpha whose end is cut obliquely
// along the edge of the merge belt (waterfall edge parallel to x), so the transfer is the same
// across the whole infeed width. Since this end is not rectangular, the infeed is a static mesh
// (parallelogram prism, written to the output folder at startup) that collides as a convex
// polyhedron (MeshShape::convex), which gives full contact patches like a box. No side guides.
//
// Parcel (HS008): rigid cuboid 480 x 330 x 300 mm, 10.75 kg, uniform mass (I_zz = 0.304 kg m^2),
// friction coefficient 0.32 (config friction.default_mu), restitution 0, short side leading, start
// with the infeed velocity and zero angular velocity.
//
// Scene parameters (config keys scene.<name>): waterfall_height [m] (default 0.01),
// merge_speed [m/s] (2.2), angle_deg (30), start_distance [m] (0.2, distance of the parcel's
// leading corner from the waterfall edge), infeed_length [m] (2.0, length of the infeed centerline
// from the waterfall edge), infeed_width [m] (0.8).
class ConveyorMergeScene : public SceneBase {
   public:
    const char* sceneName() const override { return "conveyor_merge"; }

    void populate(physics::PhysicsEngine& engine) override {
        const auto& cfg = engine.config();
        const real_t h = cfg.sceneParam("waterfall_height", (real_t)0.01);
        const real_t vMerge = cfg.sceneParam("merge_speed", (real_t)2.2);
        const real_t alpha = cfg.sceneParam("angle_deg", (real_t)30) * (real_t)M_PI / (real_t)180;
        const real_t startDistance = cfg.sceneParam("start_distance", (real_t)0.2);
        const real_t infeedScale = cfg.sceneParam("infeed_scale", (real_t)1.0);
        const real_t vInfeed = vMerge / std::cos(alpha) * infeedScale;
        const Vector3r d(-std::cos(alpha), -std::sin(alpha), 0);  // infeed direction of motion
        const real_t mass = 10.75;

        const real_t mergeHalfWidth = (real_t)0.65;
        const real_t beltHalfThickness = (real_t)0.05;

        // Merge belt: x in [-10, 2], y in [-0.65, 0.65], surface at z = 0, moving along -x.
        const auto merge = engine.addStaticBody(physics::CubeShape(Vector3r(6.0, mergeHalfWidth, beltHalfThickness)),
                                                physics::RigidState(Vector3r(-4.0, 0.0, -beltHalfThickness)));
        engine.setSurfaceVelocity(merge, Vector3r(-vMerge, 0, 0));

        // Parcel: long side (x) along the direction of motion, i.e. short side leading.
        const Vector3r half(0.24, 0.165, 0.15);
        const real_t yaw = std::atan2(d.y(), d.x());
        // Half extent of the footprint perpendicular to the waterfall edge (along y).
        const real_t footprintHalfY = half.x() * std::abs(std::sin(yaw)) + half.y() * std::abs(std::cos(yaw));
        const Vector3r parcelStart(0.0, mergeHalfWidth + startDistance + footprintHalfY, h + half.z());

        // Infeed belt: strip along -d through the parcel's start, width w, surface at z = h, end cut
        // along the waterfall edge y = 0.65.
        const real_t w = cfg.sceneParam("infeed_width", (real_t)0.8);
        const real_t length = cfg.sceneParam("infeed_length", (real_t)2.0);
        const std::string infeedMesh = writeInfeedMesh(cfg.output_folder + "/" + cfg.output_filename_prefix + "_infeed.obj", parcelStart, -d, mergeHalfWidth, w, length, h,
                                                       2 * beltHalfThickness);
        physics::MeshShape infeedShape(infeedMesh);
        infeedShape.convex = true;
        const auto infeed = engine.addRigidBody(infeedShape, physics::RigidState(Vector3r::Zero()), physics::RigidProps{});  // no mass: static
        engine.setSurfaceVelocity(infeed, vInfeed * d);

        physics::RigidState state;
        state.position = parcelStart;
        state.setOrientation(Quaternion4r(Eigen::AngleAxis<real_t>(yaw, Vector3r::UnitZ())));
        state.linearVelocity = vInfeed * d;
        state.angularVelocity = Vector3r::Zero();
        m_parcel = engine.addRigidBody(physics::CubeShape(half), state, physics::RigidProps(mass));
        engine.track(m_parcel, "parcel");

        auto Idiag = engine.ecs().get<C_InertiaDiag>(m_parcel).I;
        std::cout << "Inertia diag: " << Idiag.transpose() << std::endl;
    }

    entt::entity parcel() const { return m_parcel; }

   private:
    // Writes the infeed as a parallelogram prism (OBJ, world coordinates) and returns its path. The
    // centerline runs through `through` along the upstream direction `a` (unit, in the xy-plane);
    // the downstream end is cut along the line y = yEdge, the upstream end perpendicular to a at
    // the distance `length` from the edge (measured along the centerline).
    static std::string writeInfeedMesh(const std::string& path, const Vector3r& through, const Vector3r& a, real_t yEdge, real_t width, real_t length, real_t zTop, real_t thickness) {
        const Vector3r n(-a.y(), a.x(), 0);                                              // in-plane normal of the centerline
        const Vector3r p0 = through - ((through.y() - yEdge) / a.y()) * a;               // centerline on the edge
        // Outline: along the edge, up one lateral line, across the upstream end, back down.
        const real_t hw = (real_t)0.5 * width;
        const auto onEdge = [&](real_t s) { return Vector3r(p0 + s * n - (s * n.y() / a.y()) * a); };  // lateral line meets y = yEdge
        std::array<Vector3r, 4> corners = {onEdge(-hw), onEdge(hw), Vector3r(p0 + length * a + hw * n), Vector3r(p0 + length * a - hw * n)};
        // Make the outline counterclockwise in the top view (outward-oriented faces below).
        const Vector3r e0 = corners[1] - corners[0], e1 = corners[2] - corners[1];
        if (e0.x() * e1.y() - e0.y() * e1.x() < 0) std::swap(corners[1], corners[3]);

        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
        std::ofstream f(path);
        f.precision(12);
        for (const auto& c : corners) f << "v " << c.x() << " " << c.y() << " " << zTop - thickness << "\n";  // 1..4 bottom
        for (const auto& c : corners) f << "v " << c.x() << " " << c.y() << " " << zTop << "\n";              // 5..8 top
        // Outward-oriented triangles (counterclockwise seen from outside).
        f << "f 5 6 7\nf 5 7 8\n";  // top
        f << "f 1 3 2\nf 1 4 3\n";  // bottom
        for (int k = 0; k < 4; ++k) {  // sides
            const int b0 = 1 + k, b1 = 1 + (k + 1) % 4, t0 = b0 + 4, t1 = b1 + 4;
            f << "f " << b0 << " " << b1 << " " << t1 << "\nf " << b0 << " " << t1 << " " << t0 << "\n";
        }
        return path;
    }

    entt::entity m_parcel{entt::null};
};
