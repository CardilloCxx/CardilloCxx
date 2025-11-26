#pragma once

#include <entt/entt.hpp>
#include <optional>
#include <cmath>
#include <limits>
#include "../misc/types.hpp"
#include "physics_system.hpp"
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

// Joint-frame description shared by all constraints it always attached to Body A.
struct JointProperties {
    JointProperties() = default;


    explicit JointProperties(entt::registry& reg, entt::entity A, entt::entity B) 
        : JointProperties(reg, A, B, std::nullopt, 
            reg.get<cardillo::PhysicsSystem::C_Position3>(A).value, 
            reg.get<cardillo::PhysicsSystem::C_Position3>(B).value, 
            Matrix33r::Identity()) {}

    explicit JointProperties(entt::registry& reg,
                            entt::entity A,
                            entt::entity B,
                            std::optional<entt::entity> ref,
                            const Vector3r& r_refJ1,
                            const Vector3r& r_refJ2,
                            const Matrix33r& A_refJ) 
        :entityA(A), entityB(B)
        {
            Vector3r  r_ORef;
            Matrix33r A_Ref;

            if (ref.has_value()) {
                r_ORef = reg.get<cardillo::PhysicsSystem::C_Position3>(*ref).value;
                A_Ref  = reg.get<cardillo::PhysicsSystem::C_Orientation>(*ref).value.toRotationMatrix();
            } else {
                r_ORef = Vector3r::Zero();
                A_Ref  = Matrix33r::Identity();
            }
           
            Vector3r r_OJ1_world = r_ORef + A_Ref * r_refJ1;
            Vector3r r_OJ2_world = r_ORef + A_Ref * r_refJ2;

            const auto& r_OA = reg.get<cardillo::PhysicsSystem::C_Position3>(A).value;
            const auto& A_A = reg.get<cardillo::PhysicsSystem::C_Orientation>(A).value.toRotationMatrix();
            K1_r_S1J1 = A_A.transpose() * (r_OJ1_world - r_OA);

            const auto& r_OB = reg.get<cardillo::PhysicsSystem::C_Position3>(B).value;
            const auto& A_B = reg.get<cardillo::PhysicsSystem::C_Orientation>(B).value.toRotationMatrix();
            K2_r_S2J2 = A_B.transpose() * (r_OJ2_world - r_OB);

            const auto& A_IJ = A_Ref * A_refJ;
            A_K1J = A_A.transpose() * A_IJ;  // J -> Ref -> I -> A ==== J -> A
            A_K2J = A_B.transpose() * A_IJ;  // J -> Ref -> I -> B ==== J -> B

            g = A_IJ.transpose() * (r_OJ2_world - r_OJ1_world);
        } 

    // Joint position and orientation in the chosen reference frame
    Vector3r K1_r_S1J1{Vector3r::Zero()}; // Joint pos 1 in body A frame from S1
    Vector3r K2_r_S2J2{Vector3r::Zero()}; // Joint pos 2 in body B frame from S2
    Vector3r g{Vector3r::Zero()};          // Initial Joint offset in joint frame
    Matrix33r A_K1J{Matrix33r::Identity()};  // Joint orientation that maps J -> A
    Matrix33r A_K2J{Matrix33r::Identity()};

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
                                  const JointProperties& jointProps,
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
    void buildJointJacobian(const WorldAttachments& wa,
                            MatrixXXr& WgA,
                            MatrixXXr& WgB) const;
};

// Fully rigid joint: all 6 DOFs locked (very stiff translational and rotational)
class RigidConstraint : public TranslationRotationConstraint {
public:
    RigidConstraint(entt::registry& reg,
                    entt::entity a,
                    entt::entity b)
        : TranslationRotationConstraint(reg,
                                        a,
                                        b,
                                        JointProperties(reg, a, b),
                                        Vector3r::Constant(std::numeric_limits<real_t>::infinity()), // K_trans
                                        Vector3r::Zero(),                                            // D_trans
                                        Vector3r::Constant(std::numeric_limits<real_t>::infinity()), // K_rot
                                        Vector3r::Zero()) {}                                         // D_rot
};

// Purely translational joint: configure translational stiffness; rotations free by default
class TranslationalConstraint : public TranslationRotationConstraint {
public:
    TranslationalConstraint(entt::registry& reg,
                            entt::entity a,
                            entt::entity b,
                            const JointProperties& jointProps,
                            const Vector3r& K_trans,
                            const Vector3r& D_trans = Vector3r::Zero())
        : TranslationRotationConstraint(reg, a, b, jointProps, K_trans, D_trans,
                                        /*K_rot*/ Vector3r::Zero(),
                                        /*D_rot*/ Vector3r::Zero()) {}
};

// Hinge joint: 3 translations locked, 2 rotations locked, 1 rotation free
class HingeConstraint : public TranslationRotationConstraint {
public:
    HingeConstraint(entt::registry& reg,
                    entt::entity a,
                    entt::entity b,
                    const JointProperties& jointProps,
                    const Vector3r& hingeAxisLocalA = Vector3r(0,0,1),
                    real_t K_axis = 0,
                    real_t D_axis = (real_t)0,
                    const Vector3r& K_trans = Vector3r::Constant(std::numeric_limits<real_t>::infinity()),
                    const Vector3r& D_trans = Vector3r::Zero())
        : TranslationRotationConstraint(reg, a, b, jointProps,
                                        K_trans, D_trans,
                                        /*K_rot*/ Vector3r(K_axis, std::numeric_limits<real_t>::infinity(), std::numeric_limits<real_t>::infinity()),
                                        /*D_rot*/ Vector3r(D_axis, 0, 0))
        , m_hingeAxisA(hingeAxisLocalA)
        , m_K_axis(K_axis)
        , m_D_axis(D_axis) {}

private:
    Vector3r m_hingeAxisA{Vector3r(0,0,1)};
    real_t   m_K_axis{0};
    real_t   m_D_axis{0};
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
                   const Vector3r& Df = Vector3r::Zero());

    void setKg(const Vector3r& Kg) { m_Ke = Kg; }
    void setKf(const Vector3r& Kf) { m_Kf = Kf; }
    void setDg(const Vector3r& Dg) { m_De = Dg; }
    void setDf(const Vector3r& Df) { m_Df = Df; }

    // Fill W and compliance vectors for 6 beam rows
    ConstraintResult getConstraint() const override;

private:
    Vector3r m_Ke{Vector3r::Zero()};
    Vector3r m_Kf{Vector3r::Zero()};
    Vector3r m_De{Vector3r::Zero()};
    Vector3r m_Df{Vector3r::Zero()};
    Vector3r m_delta0{Vector3r::Zero()};
    Vector3r m_kappa0{Vector3r::Zero()};
};

} // namespace physics
} // namespace cardillo
