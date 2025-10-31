#pragma once

#include <vector>

// Project types
#include "../misc/types.hpp"
#include "physics_system.hpp"
#include <entt/entt.hpp>

namespace cardillo {
namespace physics {

// Convenient fixed-size aliases used by this header
using Vector6r = cardillo::Vectorr<6>;
using Matrix3r = cardillo::Matrix33r;
using Matrix6r = cardillo::Matrixr<6,6>;
using Matrix6x12r = cardillo::Matrixr<6,12>;

/// A single 6-DOF generalized spring constraint between two rigid bodies.
struct SpringConstraint
{
    // Bodies involved (-1 means static / world)
    entt::entity bodyA = entt::null;
    entt::entity bodyB = entt::null;

    // Local attachment points in body frames
    Vector3r rA_local = Vector3r::Zero();
    Vector3r rB_local = Vector3r::Zero();

    // Rest configuration. We store a small-rotation parametrization
    // for the rotational part (e.g. exponential map) and a translation.
    Vector3r rest_translation = Vector3r::Zero();
    Vector3r rest_rotation = Vector3r::Zero();

    // Stiffness and damping (6x6 each). For convenience default to zero.
    // World-frame stiffness/damping.
    Matrix6r K = Matrix6r::Zero();
    Matrix6r D = Matrix6r::Zero();

    // Local (body-A) stiffness/damping storage. When the user adds a spring
    // via the add*Spring helpers the axis is recorded in body-A local frame
    // and written into these local K/D matrices. These are transformed each
    // update into world-frame K/D using the body-A orientation RA.
    Matrix6r K_local = Matrix6r::Zero();
    Matrix6r D_local = Matrix6r::Zero();

    // Cached world-space transforms and attachment positions
    Matrix3r RA = Matrix3r::Identity();
    Matrix3r RB = Matrix3r::Identity();
    Vector3r xA = Vector3r::Zero();
    Vector3r xB = Vector3r::Zero();

    // Derived deformation vectors (6x1): [pos_error; rot_error]
    Vector6r g = Vector6r::Zero();
    Vector6r gdot = Vector6r::Zero();

    // Each constraint holds a pointer to the ECS registry provided at construction
    // The constructor reads the current
    // poses from the registry to initialize rest_translation/rest_rotation.
    SpringConstraint(entt::registry& reg,
                                entt::entity a,
                                entt::entity b,
                                const Vector3r& rA_loc = Vector3r::Zero(),
                                const Vector3r& rB_loc = Vector3r::Zero());

    // Pointer to the registry (non-owning). Must remain valid for the
    // lifetime of this constraint instance.
    entt::registry* registry = nullptr;
    
    // Add a translational spring along a given world-space axis with stiffness k and damping d
    // This axis will be rotated with BodyA as it move, so its stored in As local space. 
    // Special cases where k = infinity can be used to create hard constraints along that axis.
    // When k = 0, no spring is added along that axis.
    void addTranslationalSpring(const Vector3r& axis_world, real_t k, real_t d);

    // Add a rotational spring around a given world-space axis with stiffness k and damping d
    // This axis will be rotated with BodyA as it moves, so its stored in As local space.
    void addRotationalSpring(const Vector3r& axis_world, real_t k, real_t d);

    // Helpers for special cases:
    void fixTranslation();
    void fixRotation();

    // Compute world-space attachment transforms and positions from body states.
    // This updates RA, RB, xA, xB and g/gdot (position/rotation error and rate).
    // It does NOT modify global assembler structures; it only computes local
    // cached quantities per-constraint.
    void recomputeDeformations();

    // Compute Jacobian blocks for this constraint into the provided matrices.
    // Wg and Wgamma are 6x12 matrices (6 rows for constraint DOFs, 12 columns for
    // the two bodies' 6 DOFs each). They are zero-initialized on input.
    void computeWBlocks(Matrix6r & WgA_out, Matrix6r & WgB_out,
                        Matrix6r & WgammaA_out, Matrix6r & WgammaB_out) const;

};

/// Manager that stores a list of spring constraints and updates them.
class ForceInteractionManager
{
public:
    ForceInteractionManager() = default;

    // Add a constraint. Returns index of the newly added constraint.
    size_t addConstraint(const SpringConstraint& c);

    // Remove a constraint by index. If index invalid, no-op.
    void removeConstraint(size_t idx);

    // Accessors
    const std::vector<SpringConstraint>& constraints() const { return m_constraints; }
    std::vector<SpringConstraint>& constraints() { return m_constraints; }

    // Update all constraints given current rigid-body states (array-like).
    // The states container must be indexable by body index: states[b] -> RigidBodyState.
    void updateAll()
    {
        for (auto &c : m_constraints) {
            c.recomputeDeformations();
        }
    }

private:
    std::vector<SpringConstraint> m_constraints;
};

} // namespace physics
} // namespace cardillo
