#pragma once

#include <entt/entt.hpp>
#include <optional>
#include <cmath>
#include <limits>
#include "../misc/types.hpp"
#include "physics_system.hpp"

namespace cardillo {
namespace physics {

// Aliases for compact Jacobian rows (1 x 6 per body)
using Row6r = Matrixr<1,6>;

// Multi-row constraint result with dense row-storage for cleaner assembly
struct ConstraintResult {
    entt::entity a{entt::null};
    entt::entity b{entt::null};
    // Each W* is N x 6 (rows = scalar constraint rows)
    MatrixXXr WgA;     
    MatrixXXr WgB;     
    MatrixXXr WgammaA; 
    MatrixXXr WgammaB; 
    VectorXr Crows;    // size N (compliances per spring row)
    VectorXr Arows;    // size N (compliances per damper row)
};

// Base pattern: stores registry, entities and attachment points and shared helpers
class ConstraintPattern {
public:
    ConstraintPattern(entt::registry& reg,
                      entt::entity a,
                      entt::entity b,
                      const Vector3r& rA_local = Vector3r::Zero(),
                      const Vector3r& rB_local = Vector3r::Zero())
        : m_reg(&reg), m_a(a), m_b(b), m_rA_local(rA_local), m_rB_local(rB_local) {}

    virtual ~ConstraintPattern() = default;

    entt::entity entityA() const { return m_a; }
    entt::entity entityB() const { return m_b; }

    // Compute the constraint data ready for assembly
    virtual ConstraintResult getConstraint() const = 0;

protected:
    // Helper to fetch world transforms and attachment points
    struct WorldAttachments {
        Vector3r xA{Vector3r::Zero()};
        Vector3r xB{Vector3r::Zero()};
        Matrix33r RA{Matrix33r::Identity()};
        Matrix33r RB{Matrix33r::Identity()};
        Quaternion4r qA{Quaternion4r::Identity()};
        Quaternion4r qB{Quaternion4r::Identity()};
    };

    WorldAttachments computeAttachments_() const {
        WorldAttachments wa;
        if (m_reg) {
            Vector3r pA = Vector3r::Zero();
            Vector3r pB = Vector3r::Zero();
            Quaternion4r qA = Quaternion4r::Identity();
            Quaternion4r qB = Quaternion4r::Identity();
            if (m_a != entt::null) {
                if (m_reg->all_of<PhysicsSystem::C_Position3>(m_a)) pA = m_reg->get<PhysicsSystem::C_Position3>(m_a).value;
                if (m_reg->all_of<PhysicsSystem::C_Orientation>(m_a)) qA = m_reg->get<PhysicsSystem::C_Orientation>(m_a).q;
            }
            if (m_b != entt::null) {
                if (m_reg->all_of<PhysicsSystem::C_Position3>(m_b)) pB = m_reg->get<PhysicsSystem::C_Position3>(m_b).value;
                if (m_reg->all_of<PhysicsSystem::C_Orientation>(m_b)) qB = m_reg->get<PhysicsSystem::C_Orientation>(m_b).q;
            }
            wa.qA = qA;
            wa.qB = qB;
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
    Vector3r m_rA_local{Vector3r::Zero()};
    Vector3r m_rB_local{Vector3r::Zero()};
};

// Linear distance-based spring along current A->B direction.
class LinearDistanceConstraint : public ConstraintPattern {
public:
    LinearDistanceConstraint(entt::registry& reg,
                             entt::entity a,
                             entt::entity b,
                             const Vector3r& rA_local = Vector3r::Zero(),
                             const Vector3r& rB_local = Vector3r::Zero(),
                             real_t stiffness = (real_t)0,
                             real_t damping = (real_t)0)
        : ConstraintPattern(reg, a, b, rA_local, rB_local), m_k(stiffness), m_d(damping) {}

    void setStiffness(real_t k) { m_k = k; }
    void setDamping(real_t d) { m_d = d; }

    ConstraintResult getConstraint() const override {
        ConstraintResult out; out.a = m_a; out.b = m_b;

        const auto wa = computeAttachments_();
        const Vector3r r = wa.xB - wa.xA;
        const real_t len = r.norm();
        Vector3r n = (len > (real_t)0) ? (r / len) : Vector3r((real_t)1, (real_t)0, (real_t)0);

        const Vector3r rAw = wa.RA * m_rA_local;
        const Vector3r rBw = wa.RB * m_rB_local;
        // Prepare single-row matrices (1 x 6)
        out.WgA = MatrixXXr::Zero(1, 6);
        out.WgB = MatrixXXr::Zero(1, 6);
        out.WgammaA = MatrixXXr::Zero(1, 6);
        out.WgammaB = MatrixXXr::Zero(1, 6);
        // d g / d vA = -n^T, d g / d wA = (rAw x n)^T
        out.WgA.block(0, 0, 1, 3) = -n.transpose();
        out.WgA.block(0, 3, 1, 3) = (rAw.cross(n)).transpose();
        // d g / d vB = +n^T, d g / d wB = -(rBw x n)^T
        out.WgB.block(0, 0, 1, 3) =  n.transpose();
        out.WgB.block(0, 3, 1, 3) = (-rBw.cross(n)).transpose();

        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Compliance (size 1)
        out.Crows = VectorXr::Zero(1);
        out.Arows = VectorXr::Zero(1);
        const real_t eps = (real_t)1e-12;
        // Allow infinite stiffness/damping (C=0/A=0) when user passes +inf
        if (std::isfinite(m_k)) out.Crows[0] = (m_k > eps) ? (real_t)1 / m_k : (real_t)1 / eps;
        else out.Crows[0] = (real_t)0;

        if (std::isfinite(m_d)) out.Arows[0] = (m_d > eps) ? (real_t)1 / m_d : (real_t)1 / eps;
        else out.Arows[0] = (real_t)0;

        return out;
    }

private:
    real_t m_k{(real_t)0};
    real_t m_d{(real_t)0};
};

// Beam constraint (6 scalar rows): stretch/shear (x,y,z) and torsion/bend (x,y,z)
// Ke: [E1, E2, E3] stiffnesses for generalized stretch/shear
// Kf: [F1, F2, F3] stiffnesses for torsion/bending
class BeamConstraint : public ConstraintPattern {
public:
    BeamConstraint(entt::registry& reg,
                   entt::entity a,
                   entt::entity b,
                   const Vector3r& rA_local,
                   const Vector3r& rB_local,
                   const Vector3r& Kg,
                   const Vector3r& Kf,
                   const Vector3r& Dg = Vector3r::Zero(),
                   const Vector3r& Df = Vector3r::Zero())
        : ConstraintPattern(reg, a, b, rA_local, rB_local), m_Ke(Kg), m_Kf(Kf), m_De(Dg), m_Df(Df) 
        {
            const auto wa = computeAttachments_();
            m_delta0 = wa.xB - wa.xA;
        }

    void setKg(const Vector3r& Kg) { m_Ke = Kg; }
    void setKf(const Vector3r& Kf) { m_Kf = Kf; }
    void setDg(const Vector3r& Dg) { m_De = Dg; }
    void setDf(const Vector3r& Df) { m_Df = Df; }

    // Fill W and compliance vectors for 6 beam rows
    ConstraintResult getConstraint() const override {
        ConstraintResult out; out.a = m_a; out.b = m_b;
        const auto wa = computeAttachments_();

        // Mid-orientation and strains (should we use Eigens slerp instead?)
        // Quaternion4r qMid = (real_t)0.5 * (wa.qA + wa.qB);  (Addition not directly defined for quaternions)
        // Quaternion4r qMid( 
        //     (real_t)0.5 * (wa.qA.w() + wa.qB.w()),
        //     (real_t)0.5 * (wa.qA.x() + wa.qB.x()),
        //     (real_t)0.5 * (wa.qA.y() + wa.qB.y()),
        //     (real_t)0.5 * (wa.qA.z() + wa.qB.z())
        // );

        const Quaternion4r qMid = wa.qB.slerp((real_t)0.5, wa.qA);
        const Matrix33r A_mid = qMid.toRotationMatrix();
        const Vector3r gamma = A_mid.transpose() * (wa.xB - wa.xA - m_delta0);   // axial + shear strain
        const Matrix33r gamma_skew = skew_from_vector(gamma);

        // What about Ke^-1 n_i + gamma_i 0 = e_x ?

        // Differential orientation / curvature (should we use proper quaternion diff, like wa.qB.conjugate() * wa.qA?)
        const Vector3r Qd_q(
            (wa.qB.x() - wa.qA.x()),
            (wa.qB.y() - wa.qA.y()),
            (wa.qB.z() - wa.qA.z())
        );
        const real_t   Qd_w = (wa.qB.w() - wa.qA.w());

        // Quaternion4r Q_rel = wa.qB * wa.qA.conjugate();
        // Q_rel = Quaternion4r::Identity();                   Doing this does not change anything in the solution?
        // const Vector3r Qd_q( Q_rel.x(), Q_rel.y(), Q_rel.z() );
        // const real_t   Qd_w = Q_rel.w();

        const Vector3r Qmid_q(qMid.x(), qMid.y(), qMid.z());
        const real_t   Qmid_w = qMid.w();
        const real_t   factor = (real_t)2.0 / qMid.squaredNorm();
        const Vector3r kappa = factor * (Qmid_w * Qd_q - Qmid_q.cross(Qd_q) - Qmid_q * Qd_w);  // curvature
        const Matrix33r kappa_skew = skew_from_vector(kappa);

        out.WgB = MatrixXXr::Zero(6, 6);
        out.WgB.block(0, 0, 3, 3) =  A_mid;
        out.WgB.block(3, 0, 3, 3) =  (real_t) -0.5 * gamma_skew;
        out.WgB.block(3, 3, 3, 3) =  Matrix33r::Identity() - (real_t) 0.5 * kappa_skew;  // l_0 jeweils rausgekürzt

        out.WgA = -out.WgB;
        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Row-wise compliances (spring and damper)
        out.Crows = VectorXr::Zero(6);
        out.Arows = VectorXr::Zero(6);

        const real_t eps = (real_t)1e-12;
        auto asCompliance = [&](real_t k)->real_t{
            // Setting C, A to infinity prevents the corresponding constraints from being inserted into the system
            if (std::isfinite(k)) return (k > eps) ? (real_t)1 / k : std::numeric_limits<real_t>::infinity(); else return (real_t)0;
        };

        out.Crows[0] = asCompliance(m_Ke.x()); out.Crows[1] = asCompliance(m_Ke.y()); out.Crows[2] = asCompliance(m_Ke.z());
        out.Crows[3] = asCompliance(m_Kf.x()); out.Crows[4] = asCompliance(m_Kf.y()); out.Crows[5] = asCompliance(m_Kf.z());

        out.Arows[0] = asCompliance(m_De.x()); out.Arows[1] = asCompliance(m_De.y()); out.Arows[2] = asCompliance(m_De.z());
        out.Arows[3] = asCompliance(m_Df.x()); out.Arows[4] = asCompliance(m_Df.y()); out.Arows[5] = asCompliance(m_Df.z());

        return out;
    }

private:
    Vector3r m_Ke{Vector3r::Zero()};
    Vector3r m_Kf{Vector3r::Zero()};
    Vector3r m_De{Vector3r::Zero()};
    Vector3r m_Df{Vector3r::Zero()};
    Vector3r m_delta0{Vector3r::Zero()};
};

} // namespace physics
} // namespace cardillo
