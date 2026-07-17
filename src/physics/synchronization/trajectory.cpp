#include "trajectory.hpp"

#include <algorithm>
#include <cmath>

namespace cardillo {
namespace physics {
namespace {

constexpr real_t kEps = (real_t)1e-12;

TrajectoryTwist buildTwistFromEcs(const entt::registry& reg, entt::entity e) {
    Vector3r v = Vector3r::Zero();
    Vector3r w = Vector3r::Zero();
    if (reg.any_of<C_LinearVelocity3>(e)) v = reg.get<C_LinearVelocity3>(e).value;
    if (reg.any_of<C_AngularVelocity3>(e)) w = reg.get<C_AngularVelocity3>(e).value;
    return {v, w};
}

void applyPose(World& world, entt::entity e, const TrajectoryPose& pose, bool hasOrientation) {
    world.setPosition(e, pose.first);
    if (hasOrientation) {
        Quaternion4r q = pose.second;
        if (!q.coeffs().allFinite() || q.coeffs().squaredNorm() <= kEps) {
            q = Quaternion4r::Identity();
        } else {
            q.normalize();
        }
        world.setOrientation(e, q);
    }
}

void applyTwist(World& world, entt::entity e, const TrajectoryTwist& twist, bool hasAngularVelocity) {
    world.setLinearVelocity(e, twist.first);
    if (hasAngularVelocity) {
        world.setAngularVelocity(e, twist.second);
    }
}

TrajectoryTwist differentiatePose(const TrajectoryPose& prev, const TrajectoryPose& curr, real_t dt, bool hasAngularVelocity) {
    TrajectoryTwist out;
    out.first = (curr.first - prev.first) / dt;
    out.second = Vector3r::Zero();

    if (!hasAngularVelocity) return out;

    Quaternion4r q0 = prev.second;
    Quaternion4r q1 = curr.second;
    if (!q0.coeffs().allFinite() || q0.coeffs().squaredNorm() <= kEps) q0 = Quaternion4r::Identity();
    if (!q1.coeffs().allFinite() || q1.coeffs().squaredNorm() <= kEps) q1 = Quaternion4r::Identity();
    q0.normalize();
    q1 = MathHelper::alignQuaternionTo(q1.normalized(), q0);

    // calculate relative quaternion in body-fixed basis
    Quaternion4r dq = q0.conjugate() * q1;
    if (!dq.coeffs().allFinite() || dq.coeffs().squaredNorm() <= kEps) return out;
    dq.normalize();
    
    // convert to Angle-Axis representation
    AngleAxis3r angle_axis(dq);
    
    // angular velocity: scale the axis by the angular speed (angle / dt)
    out.second = angle_axis.axis() * (angle_axis.angle() / dt);
    return out;
}

}  // namespace

void Trajectory::update(World& world, real_t dt) {
    auto& reg = world.ecs();

    auto view = reg.view<C_StaticTrajectory>();
    for (auto [e, traj] : view.each()) {
        if (!reg.valid(e)) continue;

        const bool hasPos = traj.positionFunc.has_value();
        const bool hasVel = traj.velocityFunc.has_value();
        if (!hasPos && !hasVel) continue;

        const bool hasOrientation = reg.any_of<C_Orientation>(e);
        const bool hasAngularVelocity = reg.any_of<C_AngularVelocity3>(e);

        const real_t t = traj.elapsed;
        const TrajectoryTwist twist_current = buildTwistFromEcs(reg, e);

        std::optional<TrajectoryPose> pose_cmd;
        std::optional<TrajectoryTwist> twist_cmd;

        if (hasPos) {
            pose_cmd = (*traj.positionFunc)(t);
        }
        if (hasVel) {
            twist_cmd = (*traj.velocityFunc)(t);
        }

        if (hasPos && !hasVel) {
            if (dt > kEps) {
                const TrajectoryPose pose_next = (*traj.positionFunc)(t + dt);
                twist_cmd = differentiatePose(*pose_cmd, pose_next, dt, hasAngularVelocity);
            } else {
                twist_cmd = twist_current;
            }
        }

        if (hasPos && pose_cmd.has_value()) {
            applyPose(world, e, *pose_cmd, hasOrientation);
            traj.previousPosition = *pose_cmd;
        } else {
            // Velocity-authoritative trajectories should be integrated by the global integrator.
            traj.previousPosition.reset();
        }

        if (twist_cmd.has_value()) {
            applyTwist(world, e, *twist_cmd, hasAngularVelocity);
        }

        traj.initialized = true;
        traj.elapsed += std::max((real_t)0, dt);
    }
}

}  // namespace physics
}  // namespace cardillo
