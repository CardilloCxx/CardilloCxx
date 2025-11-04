#pragma once

#include <entt/entt.hpp>
#include <optional>
#include <cmath>
#include "../misc/types.hpp"
#include "physics_system.hpp"

namespace cardillo {
namespace physics {

// Aliases for compact Jacobian rows (1 x 6 per body)
using Row6r = cardillo::Matrixr<1,6>;
using Row1r = cardillo::Matrixr<1,1>;

struct ConstraintResult {
    // Entities constrained
    entt::entity a{entt::null};
    entt::entity b{entt::null};

    // Jacobian rows (along the constrained scalar direction)
    Row6r WgA;     // [dvA(3) | dwA(3)]
    Row6r WgB;     // [dvB(3) | dwB(3)]
    Row6r WgammaA; // viscous counterpart
    Row6r WgammaB; // viscous counterpart

    // Compliance matrices for stiffness (C = K^{-1}) and damping (A = D^{-1})
    // For scalar distance constraints these are 1x1.
    Row1r C;   // compliance (1/k)
    Row1r A;   // damping compliance (1/d)

    // Optional
    real_t restLength{(real_t)0};
    real_t currentLength{(real_t)0};
    cardillo::Vector3r dir_n = cardillo::Vector3r::UnitX(); // unit direction from A->B
};

// Base pattern: stores registry, entities and attachment points and shared helpers
class ConstraintPattern {
public:
    ConstraintPattern(entt::registry& reg,
                      entt::entity a,
                      entt::entity b,
                      const cardillo::Vector3r& rA_local = cardillo::Vector3r::Zero(),
                      const cardillo::Vector3r& rB_local = cardillo::Vector3r::Zero())
        : m_reg(&reg), m_a(a), m_b(b), m_rA_local(rA_local), m_rB_local(rB_local) {}

    virtual ~ConstraintPattern() = default;

    entt::entity entityA() const { return m_a; }
    entt::entity entityB() const { return m_b; }

    // Compute the constraint data ready for assembly
    virtual ConstraintResult getConstraint() const = 0;

protected:
    // Helper to fetch world transforms and attachment points
    struct WorldAttachments {
        cardillo::Vector3r xA{cardillo::Vector3r::Zero()};
        cardillo::Vector3r xB{cardillo::Vector3r::Zero()};
        cardillo::Matrix33r RA{cardillo::Matrix33r::Identity()};
        cardillo::Matrix33r RB{cardillo::Matrix33r::Identity()};
    };

    WorldAttachments computeAttachments_() const {
        WorldAttachments wa;
        if (m_reg) {
            cardillo::Vector3r pA = cardillo::Vector3r::Zero();
            cardillo::Vector3r pB = cardillo::Vector3r::Zero();
            cardillo::Quaternion4r qA = cardillo::Quaternion4r::Identity();
            cardillo::Quaternion4r qB = cardillo::Quaternion4r::Identity();
            if (m_a != entt::null) {
                if (m_reg->all_of<cardillo::PhysicsSystem::C_Position3>(m_a)) pA = m_reg->get<cardillo::PhysicsSystem::C_Position3>(m_a).value;
                if (m_reg->all_of<cardillo::PhysicsSystem::C_Orientation>(m_a)) qA = m_reg->get<cardillo::PhysicsSystem::C_Orientation>(m_a).q;
            }
            if (m_b != entt::null) {
                if (m_reg->all_of<cardillo::PhysicsSystem::C_Position3>(m_b)) pB = m_reg->get<cardillo::PhysicsSystem::C_Position3>(m_b).value;
                if (m_reg->all_of<cardillo::PhysicsSystem::C_Orientation>(m_b)) qB = m_reg->get<cardillo::PhysicsSystem::C_Orientation>(m_b).q;
            }
            wa.RA = qA.toRotationMatrix();
            wa.RB = qB.toRotationMatrix();
            wa.xA = pA + wa.RA * m_rA_local;
            wa.xB = pB + wa.RB * m_rB_local;
        }
        return wa;
    }

    entt::registry* m_reg{nullptr};
    entt::entity m_a{entt::null};
    entt::entity m_b{entt::null};
    cardillo::Vector3r m_rA_local{cardillo::Vector3r::Zero()};
    cardillo::Vector3r m_rB_local{cardillo::Vector3r::Zero()};
};

// Linear distance-based spring along current A->B direction.
class LinearDistanceConstraint : public ConstraintPattern {
public:
    LinearDistanceConstraint(entt::registry& reg,
                             entt::entity a,
                             entt::entity b,
                             const cardillo::Vector3r& rA_local = cardillo::Vector3r::Zero(),
                             const cardillo::Vector3r& rB_local = cardillo::Vector3r::Zero(),
                             real_t stiffness = (real_t)0,
                             real_t damping = (real_t)0)
        : ConstraintPattern(reg, a, b, rA_local, rB_local), m_k(stiffness), m_d(damping) {}

    void setStiffness(real_t k) { m_k = k; }
    void setDamping(real_t d) { m_d = d; }

    ConstraintResult getConstraint() const override {
        ConstraintResult out; out.a = m_a; out.b = m_b;

        const auto wa = computeAttachments_();
        const cardillo::Vector3r r = wa.xB - wa.xA;
        const real_t len = r.norm();
        cardillo::Vector3r n = (len > (real_t)0) ? (r / len) : cardillo::Vector3r((real_t)1, (real_t)0, (real_t)0);

        const cardillo::Vector3r rAw = wa.RA * m_rA_local;
        const cardillo::Vector3r rBw = wa.RB * m_rB_local;
        // Fill Jacobian rows
        out.WgA.setZero(); out.WgB.setZero(); out.WgammaA.setZero(); out.WgammaB.setZero();
        // d g / d vA = -n^T, d g / d wA = (rAw x n)^T
        out.WgA.block<1,3>(0,0) = -n.transpose();
        out.WgA.block<1,3>(0,3) = (rAw.cross(n)).transpose();
        // d g / d vB = +n^T, d g / d wB = -(rBw x n)^T
        out.WgB.block<1,3>(0,0) =  n.transpose();
        out.WgB.block<1,3>(0,3) = (-rBw.cross(n)).transpose();

        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Compliance
        const real_t eps = (real_t)1e-12;
        // Allow infinite stiffness/damping (C=0/A=0) when user passes +inf
        if (std::isfinite(m_k)) out.C(0,0) = (m_k > eps) ? (real_t)1 / m_k : (real_t)1 / eps;
        else out.C(0,0) = (real_t)0;

        if (std::isfinite(m_d)) out.A(0,0) = (m_d > eps) ? (real_t)1 / m_d : (real_t)1 / eps;
        else out.A(0,0) = (real_t)0;

        out.restLength = m_restLength;
        out.currentLength = len;
        out.dir_n = n;
        return out;
    }

    void setRestLength(real_t L) { m_restLength = L; }

private:
    real_t m_k{(real_t)0};
    real_t m_d{(real_t)0};
    real_t m_restLength{(real_t)0};
};

} // namespace physics
} // namespace cardillo
