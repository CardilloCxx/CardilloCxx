#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>

using namespace cardillo;

// Forced torsional oscillator: a rigid bar on a revolute joint (hinge) to the static world with a
// linear torsional spring k about the hinge axis, driven by the external moment M(t) = M0 sin(W t)
// about the same axis. The hinge axis is the body z-axis and passes through a point at distance d
// from the center of mass, so that the translational joint rows carry the centripetal force.
// Without gravity, the rotation angle phi about the hinge axis obeys
//
//   J phi'' + k phi = M0 sin(W t),  phi(0) = phi'(0) = 0,  J = I_zz + m d^2,
//
// with the closed-form solution (w0 = sqrt(k / J) != W)
//
//   phi(t) = M0 / (J (w0^2 - W^2)) (sin(W t) - (W / w0) sin(w0 t)).
//
// The parameters are chosen such that the bar rotates by more than 180 degrees, which requires
// the hinge spring to be linear in the angle (not in its sine) over several revolutions.
class TorsionalOscillatorScene : public SceneBase {
   public:
    const char* sceneName() const override { return "torsional_oscillator"; }

    // Geometry and mass of the bar (cube shape with these half extents).
    const Vector3r halfExtents{(real_t)0.2, (real_t)0.05, (real_t)0.05};
    const real_t mass{(real_t)4.0};
    // Distance of the hinge axis from the center of mass (along the body x-axis).
    const real_t offset{(real_t)0.15};
    // Natural frequency w0, forcing frequency W = ratio * w0, and the static deflection M0 / k
    // of the forcing amplitude (rad). The maximum angle is about staticAngle (1 + ratio) / (1 - ratio^2).
    const real_t omega0{(real_t)(2.0 * M_PI)};
    const real_t ratio{(real_t)0.6};
    // const real_t staticAngle{(real_t)2.5};
    const real_t staticAngle{(real_t)4.5};

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        m_ground = engine.addStaticBody(physics::SphereShape((real_t)0.01), physics::RigidState(Vector3r::Zero()));

        physics::RigidState state;
        state.position = Vector3r(offset, 0, 0);  // hinge axis through the world origin
        m_bar = engine.addRigidBody(physics::CubeShape(halfExtents), state, physics::RigidProps(mass));
        engine.disableCollisionBetween(m_ground, m_bar);

        m_J = engine.getInertiaDiag(m_bar)(2) + mass * offset * offset;
        m_k = m_J * omega0 * omega0;
        m_M0 = staticAngle * m_k;

        // Hinge about the body z-axis through the point (-offset, 0, 0) of the bar, with a linear
        // torsional spring k (finite compliance of the joint x-row, which is the hinge axis).
        engine.addHingeConstraint(m_ground, m_bar, physics::JointFrame::fromAxis(Vector3r(-offset, 0, 0), Vector3r::UnitZ(), m_bar), m_k);
    }

    // The moment of step [t, t + h] is evaluated at the midpoint of the step, where the theta-point
    // of the spring force lies for theta = 1/2. The torque is given in body components; the body
    // z-axis is the hinge axis.
    void updateScene(physics::PhysicsEngine& engine, real_t t) override { engine.applyForce(m_bar, Vector3r::Zero(), Vector3r(0, 0, moment(t + (real_t)0.5 * engine.timeStep()))); }

    real_t moment(real_t t) const { return m_M0 * std::sin(ratio * omega0 * t); }

    // Closed-form rotation angle for zero initial conditions.
    real_t closedForm(real_t t) const {
        const real_t W = ratio * omega0;
        return m_M0 / (m_J * (omega0 * omega0 - W * W)) * (std::sin(W * t) - ratio * std::sin(omega0 * t));
    }

    // Rotation angle of the bar about the world z-axis in (-pi, pi]; unwrap it over the steps.
    real_t wrappedAngle(physics::PhysicsEngine& engine) const {
        const Quaternion4r q = engine.ecs().get<C_Orientation>(m_bar).value;
        const Vector3r ex = q * Vector3r::UnitX();
        return std::atan2(ex.y(), ex.x());
    }

    entt::entity bar() const { return m_bar; }
    real_t inertia() const { return m_J; }
    real_t stiffness() const { return m_k; }

   private:
    entt::entity m_ground{entt::null};
    entt::entity m_bar{entt::null};
    real_t m_J{0};
    real_t m_k{0};
    real_t m_M0{0};
};
