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

    ConstraintResult getConstraint() const override;

private:
    Vector3r m_K{Vector3r::Zero()};
    Vector3r m_D{Vector3r::Zero()};
};

class RigidConstraint : public ConstraintPattern {
public:
    RigidConstraint(entt::registry& reg,
                                entt::entity a,
                                entt::entity b,
                                const Vector3r& rA_local = Vector3r::Zero(),
                                const Vector3r& rB_local = Vector3r::Zero())
        : ConstraintPattern(reg, a, b, rA_local, rB_local) {}
    ConstraintResult getConstraint() const override;
};

class HingeConstraint : public ConstraintPattern {
public:
    HingeConstraint(entt::registry& reg,
                                entt::entity a,
                                entt::entity b,
                                const Vector3r& rA_local = Vector3r::Zero(),                                        // Attachment point in A's local frame
                                const Vector3r& rB_local = Vector3r::Zero(),                                        // Attachment point in B's local frame
                                const Vector3r& axis_localA = Vector3r(0,0,1),                                      // Hinge axis in A's local frame
                                const real_t& K_axis = 0,                                                           // Stiffness along hinge axis > 0 -> rotational spring
                                const real_t& D_axis = (real_t)0,                                                   // Damping along hinge axis > 0 -> "friction"
                                const Vector2r& Kf_A = Vector2r::Constant(std::numeric_limits<real_t>::infinity()), // Stiffnesses for the two locked rotational DOFs
                                const Vector2r& Df_A = Vector2r::Zero(),                                            // Damping for the two locked rotational DOFs
                                const Vector3r& Ke_A = Vector3r::Constant(std::numeric_limits<real_t>::infinity()), // Stiffnesses for the three locked translational DOFs
                                 const Vector3r& De_A = Vector3r::Zero());                                           // Damping for the three locked translational DOFs

    ConstraintResult getConstraint() const override;

private:
    Matrix33r m_hingeFrame{Matrix33r::Identity()};
    Vector3r m_Ke{Vector3r::Zero()};
    Vector3r m_De{Vector3r::Zero()};
    Vector3r m_Kf{Vector3r::Zero()};
    Vector3r m_Df{Vector3r::Zero()};
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
