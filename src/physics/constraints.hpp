#pragma once

#include <entt/entt.hpp>
#include <optional>
#include <cmath>
#include <limits>
#include "../misc/types.hpp"
#include "world.hpp"
#include <iostream>

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

// Reference for a joint frame: position/orientation expressed in a frame entity
// (or inertial frame if ref is null).
struct JointFrame {
    Vector3r r_refJ{Vector3r::Zero()};
    Matrix33r A_refJ{Matrix33r::Identity()};
    std::optional<entt::entity> ref{std::nullopt};

    JointFrame() = default;
    JointFrame(entt::entity ref_)
        : ref(ref_) {}

    static JointFrame fromAxis(const Vector3r& p, const Vector3r& xAxis,
                               std::optional<entt::entity> ref = std::nullopt)
    {
        JointFrame jf;
        jf.r_refJ = p;
        jf.ref = ref;

        Vector3r yAxis = Vector3r::UnitY();
        if (std::abs(xAxis.dot(yAxis)) > (real_t)0.99) {
            yAxis = Vector3r::UnitZ();
        }
        Vector3r zAxis = xAxis.cross(yAxis).normalized();
        Vector3r yAxis_orth = zAxis.cross(xAxis).normalized();
        Matrix33r R;
        R.col(0) = xAxis.normalized();
        R.col(1) = yAxis_orth;
        R.col(2) = zAxis;
        jf.A_refJ = R;
        return jf;
    }
    JointFrame(const Vector3r& p, const Matrix33r& R = Matrix33r::Identity(),
               std::optional<entt::entity> ref_ = std::nullopt)
        : r_refJ(p), A_refJ(R), ref(ref_) {}
};

// Joint-frame description shared by all constraints; attachment data for A/B.
struct JointProperties {
    JointProperties() = default;

    explicit JointProperties(entt::registry& reg,
                             entt::entity A,
                             entt::entity B,
                             const JointFrame& jf)
        : entityA(A), entityB(B)
    {
        Vector3r  r_ORef;
        Matrix33r A_Ref;

        if (jf.ref.has_value()) {
            r_ORef = reg.get<cardillo::World::C_Position3>(*jf.ref).value;
            A_Ref  = reg.get<cardillo::World::C_Orientation>(*jf.ref).value.toRotationMatrix();
        } else {
            r_ORef = Vector3r::Zero();
            A_Ref  = Matrix33r::Identity();
        }

        Vector3r r_OJ_world = r_ORef + A_Ref * jf.r_refJ;

        const auto& r_OA = reg.get<cardillo::World::C_Position3>(A).value;
        const auto& A_A  = reg.get<cardillo::World::C_Orientation>(A).value.toRotationMatrix();
        K1_r_S1J = A_A.transpose() * (r_OJ_world - r_OA);

        const auto& r_OB = reg.get<cardillo::World::C_Position3>(B).value;
        const auto& A_B  = reg.get<cardillo::World::C_Orientation>(B).value.toRotationMatrix();
        K2_r_S2J = A_B.transpose() * (r_OJ_world - r_OB);

        const auto& A_IJ = A_Ref * jf.A_refJ;
        A_K1J = A_A.transpose() * A_IJ;  // J -> Ref -> I -> A ==== J -> A

        // Precompute local skew matrices
        K1_r_S1J_skew = skew_from_vector(K1_r_S1J);
        K2_r_S2J_skew = skew_from_vector(K2_r_S2J);
    }

    // Joint position and orientation in the chosen reference frame
    Vector3r K1_r_S1J{Vector3r::Zero()}; // Joint pos in body A frame from S1
    Vector3r K2_r_S2J{Vector3r::Zero()}; // Joint pos in body B frame from S2
    Matrix33r A_K1J{Matrix33r::Identity()};  // Joint orientation that maps J -> A
    Matrix33r K1_r_S1J_skew{Matrix33r::Zero()};
    Matrix33r K2_r_S2J_skew{Matrix33r::Zero()};

    // Compute current joint offset g in joint frame from world attachments
    Vector3r compute_g(const Vector3r& r_OA, const Vector3r& r_OB, const Matrix33r& A_IK1, const Matrix33r& A_IK2) const {
        const Vector3r r_S1_world = A_IK1 * K1_r_S1J;
        const Vector3r r_OJ1_world = r_OA + r_S1_world;
        const Vector3r r_S2_world = A_IK2 * K2_r_S2J;
        const Vector3r r_OJ2_world = r_OB + r_S2_world;
        const Matrix33r A_IJ = A_IK1 * A_K1J;
        return A_IJ.transpose() * (r_OJ2_world - r_OJ1_world);
    }

    entt::entity entityA{entt::null};
    entt::entity entityB{entt::null};
};

// Base pattern: stores registry, entities and attachment points and shared helpers
class ConstraintPattern {
public:
    ConstraintPattern(entt::registry& reg,
                      entt::entity a,
                      entt::entity b,
                      const Vector3r& rA_local = Vector3r::Zero(),
                      const Vector3r& rB_local = Vector3r::Zero());

    virtual ~ConstraintPattern() = default;

    entt::entity entityA() const { return m_a; }
    entt::entity entityB() const { return m_b; }

    // Compute the constraint data ready for assembly
    virtual ConstraintResult getConstraint() const = 0;

    bool getAttachPointsWorld(Vector3r& xA, Vector3r& xB) const;

protected:
    // Helper to fetch world transforms and attachment points
    struct WorldAttachments {
        Vector3r xA{Vector3r::Zero()};
        Vector3r xB{Vector3r::Zero()};
        Vector3r rA_world{Vector3r::Zero()};
        Vector3r rB_world{Vector3r::Zero()};
        Vector3r pA{Vector3r::Zero()};
        Vector3r pB{Vector3r::Zero()};
        Matrix33r RA{Matrix33r::Identity()};
        Matrix33r RB{Matrix33r::Identity()};
        Quaternion4r qA{Quaternion4r::Identity()};
        Quaternion4r qB{Quaternion4r::Identity()};
    };

    WorldAttachments computeAttachments_() const;
    // Shared compliance helpers to avoid duplication in patterns
    static real_t stiffnessToCompliance(real_t k);
    static void fillCompliance3(VectorXr& dst, int offset, const Vector3r& K);
    static void setCompliance1(VectorXr& dst, int idx, real_t k);

    entt::registry* m_reg{nullptr};
    entt::entity m_a{entt::null};
    entt::entity m_b{entt::null};
    Vector3r m_rA_local{Vector3r::Zero()};
    Vector3r m_rB_local{Vector3r::Zero()};
};

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

    ConstraintResult getConstraint() const override;

private:
    real_t m_k{(real_t)0};
    real_t m_d{(real_t)0};
};

// Generic 6-DOF translation+rotation constraint using the joint-frame Jacobian
// Specialised constraints (hinge, planar, spherical, translational, rigid)
// configure stiffness and damping in translation/rotation directions via K and D.
class TranslationRotationConstraint : public ConstraintPattern {
public:
    TranslationRotationConstraint(entt::registry& reg,
                                  entt::entity a,
                                  entt::entity b,
                                  const JointFrame& frame,
                                  const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                  const Vector3r& D_trans = Vector3r::Zero(),
                                  const Vector3r& K_rot   = Vector3r::Zero(),
                                  const Vector3r& D_rot   = Vector3r::Zero());

    ConstraintResult getConstraint() const override;

    // Access joint-frame description for visualization/debugging
    const JointProperties& jointProperties() const { return m_joint; }

protected:
    JointProperties m_joint;
    Vector3r m_K_trans{Vector3r::Zero()};
    Vector3r m_D_trans{Vector3r::Zero()};
    Vector3r m_K_rot{Vector3r::Zero()};
    Vector3r m_D_rot{Vector3r::Zero()};

    // Build full 6x6 Jacobians for a rigid joint using world attachments.
    void buildJointJacobian(const WorldAttachments& wa, MatrixXXr& WgA, MatrixXXr& WgB) const;
};

// Beam constraint (6 scalar rows): stretch/shear (x,y,z) and torsion/bend (x,y,z)
// Ke: [E1, E2, E3] stiffnesses for generalized stretch/shear
// Kf: [F1, F2, F3] stiffnesses for torsion/bending
class BeamConstraint : public ConstraintPattern {
public:
    BeamConstraint(entt::registry& reg,
                   entt::entity a,
                   entt::entity b,
                   const cardillo::World::BeamSpringParams& springs,
                   const cardillo::World::BeamCrossSection& section);

    // Fill W and compliance vectors for 6 beam rows
    ConstraintResult getConstraint() const override;

private:
    cardillo::World::BeamSpringParams m_springs{};
    cardillo::World::BeamCrossSection m_section{};
    Vector3r m_delta0{Vector3r::Zero()};
    Vector3r m_kappa0{Vector3r::Zero()};
    real_t   l_0{ 0 };
};

} // namespace physics
} // namespace cardillo
