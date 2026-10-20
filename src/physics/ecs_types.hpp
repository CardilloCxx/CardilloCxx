#pragma once

#include <Eigen/Dense>
#include <entt/entt.hpp>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include "../misc/types.hpp"

// Public ECS component/tag types
namespace cardillo {

struct C_Mass {
    real_t m{(real_t)0};
};
struct C_Position3 {
    Vector3r value{Vector3r::Zero()};
};
struct C_LinearVelocity3 {
    Vector3r value{Vector3r::Zero()};
};
struct C_AngularVelocity3 {
    Vector3r value{Vector3r::Zero()};
};
struct C_LinearAcceleration3 {
    Vector3r value{Vector3r::Zero()};
};
struct C_AngularAcceleration3 {
    Vector3r value{Vector3r::Zero()};
};
struct C_Orientation {
    Quaternion4r value{Quaternion4r::Identity()};
    Matrix33r rotation{Matrix33r::Identity()};  // kept in sync with `value`; always update both via setValue()/fromQuaternion()

    static C_Orientation fromQuaternion(const Quaternion4r& q) { return C_Orientation{q, q.toRotationMatrix()}; }
    void setValue(const Quaternion4r& q) {
        value = q;
        rotation = q.toRotationMatrix();
    }
};
struct C_DirectorTriad {
    Matrix33r value{Matrix33r::Identity()};
};
struct C_PhysicsObject {};
struct C_PointMassTag {};
struct C_RigidBodyTag {};
struct C_RigidBodyDirectorTag {};
struct C_Plane {
    Vector3r normal{Vector3r::UnitZ()};
    Vector3r up{Vector3r::UnitY()};
    real_t sizeX{(real_t)0};
    real_t sizeY{(real_t)0};
};
struct C_Cube {
    Vector3r center{Vector3r::Zero()};
    Vector3r halfExtents{Vector3r::Zero()};
    Quaternion4r q{Quaternion4r::Identity()};
};
struct C_Capsule {
    real_t radius{(real_t)0};
    real_t halfLength{(real_t)0};
    Matrix33r localRotation{Matrix33r::Identity()};
};
struct C_Cylinder {
    real_t radius{(real_t)0};
    real_t halfLength{(real_t)0};
};
struct C_Cone {
    real_t radius{(real_t)0};
    real_t height{(real_t)0};
};
struct C_Friction {
    real_t mu{(real_t)0};
};
struct C_Restitution {
    real_t normal{(real_t)0};
    real_t tangential{(real_t)0};
};
struct C_VisualObject {};
struct C_PointVisualTag {};
struct C_PlaneVisualTag {};
struct C_CubeVisualTag {};
struct C_CapsuleVisualTag {};
struct C_CylinderVisualTag {};
struct C_ConeVisualTag {};
struct C_Collidable {};
struct C_Radius {
    real_t r{(real_t)0};
};

struct C_SoftBodyVisualTag {};
struct C_SoftBodySurface {
    std::vector<Eigen::Vector3i> triangles;
    std::vector<entt::entity> nodes;
};

struct C_BeamElement {
    std::optional<entt::entity> prev;
    std::optional<entt::entity> next;
    real_t l0{(real_t)0};
    real_t l{(real_t)0};
};

struct C_Mesh {
    std::string path;
    Vector3r scale{1, 1, 1};
    // True if the body was created dynamic: its pose then refers to the mesh's center of mass and
    // principal axes, and the normalized mesh must be used for collision and output, also after
    // the body is made static (e.g. kinematically driven by a trajectory).
    bool normalized{false};
};
struct C_MeshVisualTag {};

struct C_Collider_Cube {
    Vector3r center{Vector3r::Zero()};
    Vector3r halfExtents{Vector3r::Zero()};
    Quaternion4r q{Quaternion4r::Identity()};
};

struct C_BeamHull {
    entt::entity endA{entt::null};
    entt::entity endB{entt::null};
    std::vector<Vector2r> polygonA;
    std::vector<Vector2r> polygonB;
    float radiusA{0.0f};
    float radiusB{0.0f};

    C_BeamHull() = default;
    C_BeamHull(const std::vector<Vector2r>& polyA, float rA, const std::vector<Vector2r>& polyB, float rB)
        : polygonA(polyA), polygonB(polyB), radiusA(rA), radiusB(rB) {}
    C_BeamHull(const std::vector<Vector2r>& poly, float r) : C_BeamHull(poly, r, poly, r) {}
};

struct C_Collider_BeamHull {
    entt::entity endA{entt::null};
    entt::entity endB{entt::null};
    std::vector<Vector2r> polygonA;
    std::vector<Vector2r> polygonB;
    float radiusA{0.0f};
    float radiusB{0.0f};

    C_Collider_BeamHull() = default;
    C_Collider_BeamHull(entt::entity a, entt::entity b, const std::vector<Vector2r>& polyA, float rA,
                        const std::vector<Vector2r>& polyB, float rB)
        : endA(a), endB(b), polygonA(polyA), polygonB(polyB), radiusA(rA), radiusB(rB) {}
    C_Collider_BeamHull(entt::entity a, entt::entity b, const std::vector<Vector2r>& poly, float r)
        : C_Collider_BeamHull(a, b, poly, r, poly, r) {}
};

struct C_Collider_Plane {
    Vector3r normal{Vector3r::UnitZ()};
    Vector3r up{Vector3r::UnitY()};
    real_t sizeX{(real_t)0};
    real_t sizeY{(real_t)0};
};

using TrajectoryPose = std::pair<Vector3r, Quaternion4r>;
using TrajectoryTwist = std::pair<Vector3r, Vector3r>;

// Kinematically driven body: the pose is prescribed as a function of time; the velocity is always
// derived from it (see Trajectory::update()), so that pose and velocity are consistent.
struct C_StaticTrajectory {
    std::function<TrajectoryPose(real_t)> positionFunc;
    real_t elapsed{(real_t)0};
};

struct C_Collider_Mesh {
    bool convex{false};  // collide as a convex polyhedron (see MeshShape::convex)
};
struct C_Collider_Sphere {};
struct C_Collider_Capsule {
    real_t radius{(real_t)0};
    real_t halfLength{(real_t)0};
    Matrix33r localRotation{Matrix33r::Identity()};
};
struct C_Collider_Cylinder {
    real_t radius{(real_t)0};
    real_t halfLength{(real_t)0};
};
struct C_Collider_Cone {
    real_t radius{(real_t)0};
    real_t height{(real_t)0};
};

struct C_ExternalForce {
    Vector3r f{Vector3r::Zero()};
};
struct C_ExternalTorque {
    Vector3r tau{Vector3r::Zero()};
};

struct C_TrackTag {
    std::string name;
};
// Prescribed velocity of the material points of a body's surface relative to the body, in body
// coordinates (e.g. a conveyor belt: the geometry is at rest, the belt surface moves). Only the
// contact laws see it: it is added to the relative contact velocity, projected onto the tangent
// plane of each contact, so it drives friction but never pushes along the contact normal.
struct C_SurfaceVelocity {
    Vector3r value{Vector3r::Zero()};
};
struct C_BodyIndex {
    int b{-1};  // index of the body in the stacked velocity vector; -1: none
};

// Optional inertia component for rigid bodies (body-frame diagonal inertia)
struct C_InertiaDiag {
    Vector3r I{Vector3r::Zero()};
};

}  // namespace cardillo
