#pragma once

#include <cmath>
#include <entt/entt.hpp>
#include <iostream>
#include <limits>
#include <optional>
#include "../../misc/types.hpp"
#include "../../rigid_body/rigid_body.hpp"
#include "../../rigid_body/transformations.hpp"
#include "../world.hpp"

namespace cardillo {
namespace physics {

// Aliases for compact Jacobian rows (1 x 6 per body)
using Row6r = Matrixr<1, 6>;

// Multi-row constraint result with dense row-storage for cleaner assembly
struct ConstraintResult {
    entt::entity a{entt::null};
    entt::entity b{entt::null};
    // Each W* is N x 6 (rows = scalar constraint rows)
    MatrixXXr WgA;
    MatrixXXr WgB;
    MatrixXXr WgammaA;
    MatrixXXr WgammaB;
    VectorXr Crows;  // size N (compliances per spring row)
    VectorXr Arows;  // size N (compliances per damper row)

    VectorXr positionError;

    std::vector<bool> c_used;
    std::vector<bool> a_used;
};

/// Describes the position and orientation of a joint in a reference frame.
/// When @p ref is empty all values are in world (inertial) frame.
/// Column 0 of @p A_refJ is the primary (hinge) axis for revolute joints.
struct JointFrame {
    /// Joint origin in the reference frame (metres).
    Vector3r r_refJ{Vector3r::Zero()};
    /// Joint orientation in the reference frame. Column 0 is the hinge/primary axis.
    Matrix33r A_refJ{Matrix33r::Identity()};
    /// Optional reference entity whose frame defines r_refJ and A_refJ. Empty = world frame.
    std::optional<entt::entity> ref{std::nullopt};

    JointFrame() = default;
    JointFrame(entt::entity ref_) : ref(ref_) {}

    static JointFrame fromAxis(const Vector3r& p, const Vector3r& xAxis, std::optional<entt::entity> ref = std::nullopt) {
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
    JointFrame(const Vector3r& p, const Matrix33r& R = Matrix33r::Identity(), std::optional<entt::entity> ref_ = std::nullopt) : r_refJ(p), A_refJ(R), ref(ref_) {}
};

// Joint-frame description shared by all constraints; attachment data for A/B.
struct JointProperties {
    JointProperties() = default;

    explicit JointProperties(entt::registry& reg, entt::entity A, entt::entity B, const JointFrame& jf) : entityA(A), entityB(B) {
        const auto inertial = RigidBody::RigidState::inertial();
        const auto refState = jf.ref.has_value() ? RigidBody::getState(reg, *jf.ref) : inertial;

        const auto stateA = RigidBody::getState(reg, A);
        const auto stateB = RigidBody::getState(reg, B);

        K1_r_S1J = transform::point(jf.r_refJ, refState, stateA);
        K2_r_S2J = transform::point(jf.r_refJ, refState, stateB);

        A_K1J = transform::rotation(jf.A_refJ, refState, stateA);

        // Precompute local skew matrices
        K1_r_S1J_skew = SkewSymmetricMatrix3r(K1_r_S1J);
        K2_r_S2J_skew = SkewSymmetricMatrix3r(K2_r_S2J);
    }

    // Joint position and orientation in the chosen reference frame
    Vector3r K1_r_S1J{Vector3r::Zero()};     // Joint pos in body A frame from S1
    Vector3r K2_r_S2J{Vector3r::Zero()};     // Joint pos in body B frame from S2
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
    ConstraintPattern(entt::registry& reg, entt::entity a, entt::entity b, index_t numRows);

    virtual ~ConstraintPattern() = default;

    entt::entity entityA() const { return m_a; }
    entt::entity entityB() const { return m_b; }

    // Compute the constraint data ready for assembly
    virtual ConstraintResult getConstraint() const = 0;

    virtual VectorXr getSource() const { return VectorXr::Zero(m_numRows); }
    virtual void setScalarVelocity(real_t v) { std::cerr << "Warning: setScalarVelocity not implemented for this constraint pattern" << std::endl; }
    virtual void setLinearVelocity(const Vector3r& v) { std::cerr << "Warning: setLinearVelocity not implemented for this constraint pattern" << std::endl; }
    virtual void setAngularVelocity(const Vector3r& w) { std::cerr << "Warning: setAngularVelocity not implemented for this constraint pattern" << std::endl; }

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

    // Local attachment points for the joint
    Vector3r m_rA_local{Vector3r::Zero()};
    Vector3r m_rB_local{Vector3r::Zero()};

    index_t m_numRows{0};
};

class LinearDistanceConstraint : public ConstraintPattern {
   public:
    LinearDistanceConstraint(entt::registry& reg, entt::entity a, entt::entity b, const Vector3r& rA_local = Vector3r::Zero(), const Vector3r& rB_local = Vector3r::Zero(),
                             real_t stiffness = std::numeric_limits<real_t>::infinity(), real_t damping = (real_t)0)
        : ConstraintPattern(reg, a, b, 1), m_k(stiffness), m_d(damping) {
        m_rA_local = rA_local;
        m_rB_local = rB_local;
        WorldAttachments wa = computeAttachments_();
        L0 = (wa.xB - wa.xA).norm();
    }

    ConstraintResult getConstraint() const override;
    VectorXr getSource() const override { return VectorXr::Constant(1, scalarVelocity); }
    void setScalarVelocity(real_t v) override { scalarVelocity = v; }

   private:
    real_t L0{(real_t)0};  // rest length
    real_t scalarVelocity = (real_t)0;
    real_t m_k{(real_t)0};
    real_t m_d{(real_t)0};
};

// Generic 6-DOF translation+rotation constraint using the joint-frame Jacobian
// Specialised constraints (hinge, planar, spherical, translational, rigid)
// configure stiffness and damping in translation/rotation directions via K and D.
class TranslationRotationConstraint : public ConstraintPattern {
   public:
    TranslationRotationConstraint(entt::registry& reg, entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                                  const Vector3r& D_trans = Vector3r::Zero(), const Vector3r& K_rot = Vector3r::Zero(), const Vector3r& D_rot = Vector3r::Zero());

    ConstraintResult getConstraint() const override;
    VectorXr getSource() const override {
        VectorXr src(6);
        src.head<3>() = m_translational_velocity;
        src.tail<3>() = m_angular_velocity;
        return src;
    }
    void setLinearVelocity(const Vector3r& v) override { m_translational_velocity = v; }
    void setAngularVelocity(const Vector3r& w) override { m_angular_velocity = w; }

    // Access joint-frame description for visualization/debugging
    const JointProperties& jointProperties() const { return m_joint; }

   protected:
    JointProperties m_joint;
    Vector3r m_K_trans{Vector3r::Zero()};
    Vector3r m_D_trans{Vector3r::Zero()};
    Vector3r m_K_rot{Vector3r::Zero()};
    Vector3r m_D_rot{Vector3r::Zero()};

    Vector3r m_translational_velocity{Vector3r::Zero()};
    Vector3r m_angular_velocity{Vector3r::Zero()};

    VectorXr getPositionError(const Vector3r& g, const ConstraintResult& res) const;
    VectorXr m_g0{VectorXr::Zero(6)};

    // Build full 6x6 Jacobians for a rigid joint using world attachments.
    void buildJointJacobian(const WorldAttachments& wa, const Vector3r& g, MatrixXXr& WgA, MatrixXXr& WgB) const;
};

// Beam constraint (6 scalar rows): stretch/shear (x,y,z) and torsion/bend (x,y,z)
// Ke: [E1, E2, E3] stiffnesses for generalized stretch/shear
// Kf: [F1, F2, F3] stiffnesses for torsion/bending
class BeamConstraint : public ConstraintPattern {
   public:
    BeamConstraint(entt::registry& reg, entt::entity a, entt::entity b, const physics::BeamSpringParams& springs, const physics::BeamCrossSection& section);

    // Fill W and compliance vectors for 6 beam rows
    ConstraintResult getConstraint() const override;

   private:
    physics::BeamSpringParams m_springs{};
    physics::BeamCrossSection m_section{};
    Vector3r m_gamma0{Vector3r::Zero()};
    Vector3r m_kappa0{Vector3r::Zero()};
    real_t l_0{0};
};

}  // namespace physics
}  // namespace cardillo
