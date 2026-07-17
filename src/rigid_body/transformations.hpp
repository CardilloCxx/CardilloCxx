#pragma once

#include <entt/entt.hpp>
#include <limits>
#include "../misc/math_helper.hpp"
#include "../misc/types.hpp"
#include "rigid_body.hpp"

/**
 * @brief Frame-transformation helpers between two rigid-body reference frames.
 *
 * Convention used throughout this file (matching cardillo::physics::RigidState): a RigidState's
 * `position`/`orientation`/`linearVelocity` are expressed in the world frame, while its
 * `angularVelocity` is expressed in its own body-fixed frame. Free function parameters named
 * `point`, `dir`, `v` are expressed in the *source* (`from`) body's local/body-fixed frame unless
 * documented otherwise, and return values are expressed in the *target* (`to`) body's local frame.
 */
namespace cardillo::transform {

namespace detail {

inline Quaternion4r sanitizeQuaternion(const Quaternion4r& q_in) {
    Quaternion4r q = q_in;
    if (!q.coeffs().allFinite()) return Quaternion4r::Identity();
    const real_t n2 = q.coeffs().squaredNorm();
    if (n2 <= std::numeric_limits<real_t>::epsilon()) return Quaternion4r::Identity();
    q.normalize();
    return q;
}

inline Vector3r bodyToInertialPoint(const RigidBody::RigidState& s, const Vector3r& r_body) {
    return s.position + s.rotation * r_body;
}
inline Vector3r inertialToBodyPoint(const RigidBody::RigidState& s, const Vector3r& r_inertial) {
    return s.rotation.transpose() * (r_inertial - s.position);
}

inline Vector3r bodyToInertialVector(const RigidBody::RigidState& s, const Vector3r& x_body) {
    return s.rotation * x_body;
}
inline Vector3r inertialToBodyVector(const RigidBody::RigidState& s, const Vector3r& x_inertial) {
    return s.rotation.transpose() * x_inertial;
}

inline Quaternion4r orientation(const Quaternion4r& q, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    const Quaternion4r q_in = sanitizeQuaternion(q);
    const Quaternion4r q_to_from = to.orientation.conjugate() * from.orientation;
    const Quaternion4r q_out = sanitizeQuaternion(q_to_from * q_in);
    return MathHelper::alignQuaternionTo(q_out, Quaternion4r::Identity());
}

inline Vector3r bodyToInertialVelocity(const RigidBody::RigidState& s, const Vector3r& r_inertial) {
    const Vector3r omega_world = bodyToInertialVector(s, s.angularVelocity);
    return s.linearVelocity + omega_world.cross(r_inertial);
}

}  // namespace detail

/// Re-expresses a point given in `from`'s local frame as a point in `to`'s local frame.
inline Vector3r point(const Vector3r& point, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    return detail::inertialToBodyPoint(to, detail::bodyToInertialPoint(from, point));
}

/// Re-expresses a free direction vector (no translation) from `from`'s local frame into `to`'s local frame.
inline Vector3r direction(const Vector3r& dir, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    return detail::inertialToBodyVector(to, detail::bodyToInertialVector(from, dir));
}

/// Re-expresses a rotation matrix defined relative to `from`'s axes as one relative to `to`'s axes.
inline Matrix33r rotation(const Matrix33r& R, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    return to.rotation.transpose() * from.rotation * R;
}

/// Re-expresses an orientation quaternion defined relative to `from`'s axes as one relative to `to`'s axes.
inline Quaternion4r orientation(const Quaternion4r& q, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    return detail::orientation(q, from, to);
}

/**
 * @brief Computes the relative linear velocity, in `to`'s body frame, of a point fixed to `from`.
 * @param v Velocity of the point relative to `from`, expressed in `from`'s body frame.
 * @param position The point's location, expressed in `from`'s body frame.
 */
inline Vector3r linearVelocity(const Vector3r& v, const Vector3r& position, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    const auto& fromFrame = from;
    const auto& toFrame = to;

    const Vector3r p_world = detail::bodyToInertialPoint(fromFrame, position);
    const Vector3r r_from_world = p_world - fromFrame.position;
    const Vector3r r_to_world = p_world - toFrame.position;

    const Vector3r v_world = detail::bodyToInertialVelocity(fromFrame, r_from_world) + detail::bodyToInertialVector(fromFrame, v);
    const Vector3r v_rel_world = v_world - detail::bodyToInertialVelocity(toFrame, r_to_world);
    return detail::inertialToBodyVector(toFrame, v_rel_world);
}

/**
 * @brief Computes the relative angular velocity, in `to`'s body frame, of a frame spinning at
 * `omega` relative to `from`.
 * @param omega Angular velocity relative to `from`, expressed in `from`'s body frame.
 */
inline Vector3r angularVelocity(const Vector3r& omega, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    const auto& fromFrame = from;
    const auto& toFrame = to;

    const Vector3r w_world = detail::bodyToInertialVector(fromFrame, fromFrame.angularVelocity + omega);
    return detail::inertialToBodyVector(toFrame, w_world) - toFrame.angularVelocity;
}

/// Re-expresses a full rigid-body state (of a body rigidly attached to `from`) relative to `to` instead.
inline RigidBody::RigidState rigidState(const RigidBody::RigidState& state, const RigidBody::RigidState& from, const RigidBody::RigidState& to) {
    RigidBody::RigidState out;
    out.position = point(state.position, from, to);
    out.orientation = orientation(state.orientation, from, to);
    out.rotation = out.orientation.toRotationMatrix();
    out.linearVelocity = linearVelocity(state.linearVelocity, state.position, from, to);

    const Vector3r omega_in_from = state.rotation * state.angularVelocity;
    const Vector3r omega_in_to = angularVelocity(omega_in_from, from, to);
    out.angularVelocity = out.rotation.transpose() * omega_in_to;
    return out;
}

}  // namespace cardillo::transform