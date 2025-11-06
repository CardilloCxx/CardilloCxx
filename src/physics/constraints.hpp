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

    bool getAttachPointsWorld(Vector3r& xA, Vector3r& xB) const {
        auto wa = computeAttachments_();
        xA = wa.xA; xB = wa.xB; return true;
    }

protected:
    // Helper to fetch world transforms and attachment points
    struct WorldAttachments {
        Vector3r xA{Vector3r::Zero()};
        Vector3r xB{Vector3r::Zero()};
        Vector3r rA_world{Vector3r::Zero()};
        Vector3r rB_world{Vector3r::Zero()};
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
            wa.rA_world = wa.RA * m_rA_local;
            wa.rB_world = wa.RB * m_rB_local;
            wa.xA = pA + wa.rA_world;
            wa.xB = pB + wa.rB_world;
        }
        return wa;
    }

    // Shared compliance helpers to avoid duplication in patterns
    static inline real_t stiffnessToCompliance(real_t k) {
        const real_t eps = (real_t)1e-12;
        if (std::isfinite(k)) return (k > eps) ? (real_t)1 / k : std::numeric_limits<real_t>::infinity();
        return (real_t)0; // +inf stiffness -> zero compliance
    }
    static inline void fillCompliance3(VectorXr& dst, int offset, const Vector3r& K) {
        if (dst.size() < offset + 3) dst.conservativeResize(offset + 3);
        dst[offset + 0] = stiffnessToCompliance(K.x());
        dst[offset + 1] = stiffnessToCompliance(K.y());
        dst[offset + 2] = stiffnessToCompliance(K.z());
    }
    static inline void setCompliance1(VectorXr& dst, int idx, real_t k) {
        if (dst.size() < idx + 1) dst.conservativeResize(idx + 1);
        dst[idx] = stiffnessToCompliance(k);
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
                             real_t stiffness = std::numeric_limits<real_t>::infinity(),
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
        // Prepare single-row matrices (6 x 1)
        out.WgA = MatrixXXr::Zero(6, 1);
        out.WgB = MatrixXXr::Zero(6, 1);
        // d g / d vA = -n^T, d g / d wA = (rAw x n)^T
        out.WgA.block(0, 0, 3, 1) = -n;
        out.WgA.block(3, 0, 3, 1) = -(m_rA_local.cross(wa.RA.transpose() * n));
        // d g / d vB = +n^T, d g / d wB = -(rBw x n)^T
        out.WgB.block(0, 0, 3, 1) =  n;
        out.WgB.block(3, 0, 3, 1) = (m_rB_local.cross(wa.RB.transpose() * n));

        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Compliance (size 1)
        out.Crows = VectorXr::Zero(1);
        out.Arows = VectorXr::Zero(1);
        setCompliance1(out.Crows, 0, m_k);
        setCompliance1(out.Arows, 0, m_d);

        return out;
    }

private:
    real_t m_k{(real_t)0};
    real_t m_d{(real_t)0};
};

class TranslationalConstraint : public ConstraintPattern {
public:
    TranslationalConstraint(entt::registry& reg,
                                   entt::entity a,
                                   entt::entity b,
                                   const Vector3r& rA_local = Vector3r::Zero(),
                                   const Vector3r& rB_local = Vector3r::Zero(),
                                   const Vector3r& K_A = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                   const Vector3r& D_A = Vector3r::Zero())
        : ConstraintPattern(reg, a, b, rA_local, rB_local), m_K(K_A), m_D(D_A) {}

    void setStiffnessA(const Vector3r& K_A) { m_K = K_A; }
    void setDampingA(const Vector3r& D_A) { m_D = D_A; }

    ConstraintResult getConstraint() const override {
        ConstraintResult out; 
        out.a = m_a; 
        out.b = m_b;

        const auto wa = computeAttachments_();

        const Matrix33r RA = wa.qA.toRotationMatrix();
        const Matrix33r RB = wa.qB.toRotationMatrix();

        const Vector3r pA = wa.xA;
        const Vector3r pB = wa.xB;

        const Vector3r rA_local = m_rA_local;
        const Vector3r rB_local = m_rB_local;

        const Matrix33r S_rA_local = skew_from_vector(rA_local);
        const Matrix33r S_rB_local = skew_from_vector(rB_local);
      
        // Initialize
        out.WgA = MatrixXXr::Zero(6, 3);
        out.WgA.block<3,3>(0,0) = -RA;
        out.WgA.block<3,3>(3,0) = -S_rA_local;

        out.WgB = MatrixXXr::Zero(6, 3);
        out.WgB.block<3,3>(0,0) =  RA;
        out.WgB.block<3,3>(3,0) =  (S_rB_local * wa.RB.transpose() * RA);

        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Compliance and damping: only 3 translational DOFs
        out.Crows = VectorXr::Zero(3);
        out.Arows = VectorXr::Zero(3);
        fillCompliance3(out.Crows, 0, m_K);
        fillCompliance3(out.Arows, 0, m_D);

        return out;
    }

private:
    Vector3r m_K{Vector3r::Zero()};
    Vector3r m_D{Vector3r::Zero()};
};

// Beam constraint (6 scalar rows): stretch/shear (x,y,z) and torsion/bend (x,y,z)
// Ke: [E1, E2, E3] stiffnesses for generalized stretch/shear
// Kf: [F1, F2, F3] stiffnesses for torsion/bending
class BeamConstraint : public ConstraintPattern {
public:
    BeamConstraint(entt::registry& reg,
                   entt::entity a,
                   entt::entity b,
                   const Vector3r& Kg,
                   const Vector3r& Kf,
                   const Vector3r& Dg = Vector3r::Zero(),
                   const Vector3r& Df = Vector3r::Zero())
        : ConstraintPattern(reg, a, b, Vector3r::Zero(), Vector3r::Zero()), m_Ke(Kg), m_Kf(Kf), m_De(Dg), m_Df(Df) 
        {
            // Check 

            const auto wa = computeAttachments_();
            const Quaternion4r qMid(0.5 * (wa.qA.coeffs() + wa.qB.coeffs()));
            const Matrix33r A_mid = qMid.normalized().toRotationMatrix();
            m_delta0 = A_mid.transpose() * (wa.xB - wa.xA);
        }

    void setKg(const Vector3r& Kg) { m_Ke = Kg; }
    void setKf(const Vector3r& Kf) { m_Kf = Kf; }
    void setDg(const Vector3r& Dg) { m_De = Dg; }
    void setDf(const Vector3r& Df) { m_Df = Df; }

    // Fill W and compliance vectors for 6 beam rows
    ConstraintResult getConstraint() const override {
        ConstraintResult out; out.a = m_a; out.b = m_b;
        const auto wa = computeAttachments_();

        // Mid-orientation and strains
        const Vector4r qMid = 0.5 * (wa.qA.coeffs() + wa.qB.coeffs());
        const Matrix33r A_mid = Quaternion4r(qMid).normalized().toRotationMatrix();
        const Vector3r gamma = A_mid.transpose() * (wa.xB - wa.xA); // - m_delta0;   // axial + shear strain
        const Matrix33r gamma_skew = skew_from_vector(gamma);

        const Vector4r dQ = wa.qB.coeffs() - wa.qA.coeffs();
        const real_t   factor = (real_t)2.0 / qMid.squaredNorm();
        const real_t   Qmid_w = qMid(3);
        const Vector3r Qmid_q = qMid.head<3>();
        const Vector3r kappa = factor * (Qmid_w * dQ.head<3>() - Qmid_q.cross(dQ.head<3>()) - Qmid_q * dQ(3));

        const Matrix33r kappa_skew = skew_from_vector(kappa);

        out.WgA = MatrixXXr::Zero(6, 6);
        out.WgA.block<3, 3>(0, 0) = -A_mid;
        out.WgA.block<3, 3>(3, 0) = (real_t) -0.5 * gamma_skew;
        out.WgA.block<3, 3>(3, 3) = -Matrix33r::Identity() - (real_t) 0.5 * kappa_skew;  // l_0 jeweils rausgekürzt

        out.WgB = MatrixXXr::Zero(6, 6);
        out.WgB.block<3, 3>(0, 0) = A_mid;
        out.WgB.block<3, 3>(3, 0) = (real_t) -0.5 * gamma_skew;
        out.WgB.block<3, 3>(3, 3) = Matrix33r::Identity() - (real_t) 0.5 * kappa_skew;  // l_0 jeweils rausgekürzt

        // out.WgA = -out.WgB;
        out.WgammaA = out.WgA;
        out.WgammaB = out.WgB;

        // Row-wise compliances (spring and damper)
        out.Crows = VectorXr::Zero(6);
        out.Arows = VectorXr::Zero(6);
        fillCompliance3(out.Crows, 0, m_Ke);
        fillCompliance3(out.Crows, 3, m_Kf);
        fillCompliance3(out.Arows, 0, m_De);
        fillCompliance3(out.Arows, 3, m_Df);

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
