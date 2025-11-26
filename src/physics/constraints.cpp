#include "constraints.hpp"

namespace cardillo { namespace physics {

// ===================== ConstraintPattern base =====================
ConstraintPattern::ConstraintPattern(entt::registry& reg,
                                     entt::entity a,
                                     entt::entity b,
                                     const Vector3r& rA_local,
                                     const Vector3r& rB_local)
    : m_reg(&reg), m_a(a), m_b(b), m_rA_local(rA_local), m_rB_local(rB_local) {}

bool ConstraintPattern::getAttachPointsWorld(Vector3r& xA, Vector3r& xB) const {
    auto wa = computeAttachments_();
    xA = wa.xA; xB = wa.xB; return true;
}

ConstraintPattern::WorldAttachments ConstraintPattern::computeAttachments_() const {
    WorldAttachments wa;
    if (m_reg) {
        Vector3r pA = Vector3r::Zero();
        Vector3r pB = Vector3r::Zero();
        Quaternion4r qA = Quaternion4r::Identity();
        Quaternion4r qB = Quaternion4r::Identity();
        if (m_a != entt::null) {
            if (m_reg->all_of<PhysicsSystem::C_Position3>(m_a)) pA = m_reg->get<PhysicsSystem::C_Position3>(m_a).value;
            if (m_reg->all_of<PhysicsSystem::C_Orientation>(m_a)) qA = m_reg->get<PhysicsSystem::C_Orientation>(m_a).value;
        }
        if (m_b != entt::null) {
            if (m_reg->all_of<PhysicsSystem::C_Position3>(m_b)) pB = m_reg->get<PhysicsSystem::C_Position3>(m_b).value;
            if (m_reg->all_of<PhysicsSystem::C_Orientation>(m_b)) qB = m_reg->get<PhysicsSystem::C_Orientation>(m_b).value;
        }
        wa.qA = qA;
        wa.qB = qB;
        wa.RA = qA.toRotationMatrix();
        wa.RB = qB.toRotationMatrix();
        wa.rA_world = wa.RA * m_rA_local;
        wa.rB_world = wa.RB * m_rB_local;
        wa.pA = pA;
        wa.pB = pB;
        wa.xA = pA + wa.rA_world;
        wa.xB = pB + wa.rB_world;
    }
    return wa;
}

real_t ConstraintPattern::stiffnessToCompliance(real_t k) {
    const real_t eps = (real_t)1e-12;
    if (std::isfinite(k)) return (k > eps) ? (real_t)1 / k : std::numeric_limits<real_t>::infinity();
    return (real_t)0; // +inf stiffness -> zero compliance
}

void ConstraintPattern::fillCompliance3(VectorXr& dst, int offset, const Vector3r& K) {
    if (dst.size() < offset + 3) dst.conservativeResize(offset + 3);
    dst[offset + 0] = stiffnessToCompliance(K.x());
    dst[offset + 1] = stiffnessToCompliance(K.y());
    dst[offset + 2] = stiffnessToCompliance(K.z());
}

void ConstraintPattern::setCompliance1(VectorXr& dst, int idx, real_t k) {
    if (dst.size() < idx + 1) dst.conservativeResize(idx + 1);
    dst[idx] = stiffnessToCompliance(k);
}

// ===================== LinearDistanceConstraint =====================
ConstraintResult LinearDistanceConstraint::getConstraint() const {
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

    out.Crows = VectorXr::Zero(1);
    out.Arows = VectorXr::Zero(1);
    setCompliance1(out.Crows, 0, m_k);
    setCompliance1(out.Arows, 0, m_d);

    return out;
}

// ===================== TranslationRotationConstraint =====================

TranslationRotationConstraint::TranslationRotationConstraint(entt::registry& reg,
                                                             entt::entity a,
                                                             entt::entity b,
                                                             const JointProperties& jointProps,
                                                             const Vector3r& K_trans,
                                                             const Vector3r& D_trans,
                                                             const Vector3r& K_rot,
                                                             const Vector3r& D_rot)
    : ConstraintPattern(reg, a, b, jointProps.K1_r_S1J1, jointProps.K2_r_S2J2)
    , m_joint(jointProps)
    , m_K_trans(K_trans)
    , m_D_trans(D_trans)
    , m_K_rot(K_rot)
    , m_D_rot(D_rot)
{}

// Helper: build full 6x6 Jacobians for a joint between A and B using
// the precomputed joint-frame geometry. This encodes a
// fully locked 6-DOF joint; specialised constraints will mask DOFs via C/A.
void TranslationRotationConstraint::buildJointJacobian(const ConstraintPattern::WorldAttachments& wa,
                                                       MatrixXXr& WgA,
                                                       MatrixXXr& WgB) const {

    WgA = MatrixXXr::Zero(6, 6);
    WgB = MatrixXXr::Zero(6, 6);

    // translations in joint frame
    const auto& A_IK1 = wa.RA;
    const auto& A_IK2 = wa.RB;
    const auto& A_K1J = m_joint.A_K1J;
    const auto& A_K2J = m_joint.A_K2J;
    const auto& r_OA  = wa.pA;
    const auto& r_OB  = wa.pB;
    const auto& r_K1_S1J1 = m_joint.K1_r_S1J1;
    const auto& r_K2_S2J2 = m_joint.K2_r_S2J2;

    const Vector3r r_S1_world = A_IK1 * m_joint.K1_r_S1J1; // K1_r_S1J1 is constant local offset
    const Vector3r r_OJ1_world = r_OA + r_S1_world;
    const Vector3r r_S2_world = A_IK2 * m_joint.K2_r_S2J2; // K2_r_S2J2 is constant local offset
    const Vector3r r_OJ2_world = r_OB + r_S2_world;
    const Matrix33r A_IJ = A_IK1 * m_joint.A_K1J; 
    const Vector3r g_current = A_IJ.transpose() * (r_OJ2_world - r_OJ1_world);

    const Matrix33r skew_g = skew_from_vector(g_current);
    const Matrix33r skew_r1 = skew_from_vector(r_K1_S1J1);
    const Matrix33r skew_r2 = skew_from_vector(r_K2_S2J2);

    // WgA.block<3,3>(0,0) = -A_IJ; 
    // WgA.block<3,3>(3,0) = -skew_r1 * A_K1J - A_K1J * skew_g;
    // WgB.block<3,3>(0,0) = A_IJ;
    // WgB.block<3,3>(3,0) = skew_r2 * A_IK2.transpose() * A_IJ;

    // translations
    //     WgA.block<3,3>(0,0) = -Matrix33r::Identity();
    //     WgA.block<3,3>(3,0) = -skew_r1 * A_IK1.transpose();
    // 
    //     WgB.block<3,3>(0,0) = Matrix33r::Identity();
    //     WgB.block<3,3>(3,0) = skew_r2 * A_IK2.transpose();

    // orientations
    // 1. x x y = z
    WgA.block<3,1>(3,3) = A_K1J.col(2);
    WgB.block<3,1>(3,3) = -A_K2J.col(2);

    // 2. y @ z = x
    WgA.block<3,1>(3,4) = A_K1J.col(0);
    WgB.block<3,1>(3,4) = -A_K2J.col(0);

    // 3. z x x = y
    WgA.block<3,1>(3,5) = A_K1J.col(1);
    WgB.block<3,1>(3,5) = -A_K2J.col(1);
}

ConstraintResult TranslationRotationConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;

    // Use standard attachment computation; joint is defined by m_rA_local/m_rB_local
    const auto wa = computeAttachments_();

    // Build full 6x6 Jacobians for a rigid joint at the current attachments.
    buildJointJacobian(wa, out.WgA, out.WgB);

    // For now, gamma rows mirror g rows
    out.WgammaA = out.WgA;
    out.WgammaB = out.WgB;

    // 6 scalar rows: 3 translational (0..2), 3 rotational (3..5)
    out.Crows = VectorXr::Zero(6);
    out.Arows = VectorXr::Zero(6);

    // Translational compliances
    fillCompliance3(out.Crows, 0, m_K_trans);
    fillCompliance3(out.Arows, 0, m_D_trans);

    // Rotational compliances
    fillCompliance3(out.Crows, 3, m_K_rot);
    fillCompliance3(out.Arows, 3, m_D_rot);

    return out;
}

// ===================== BeamConstraint =====================
BeamConstraint::BeamConstraint(entt::registry& reg,
                               entt::entity a,
                               entt::entity b,
                               const Vector3r& Kg,
                               const Vector3r& Kf,
                               const Vector3r& Dg,
                               const Vector3r& Df)
    : ConstraintPattern(reg, a, b, Vector3r::Zero(), Vector3r::Zero()), m_Ke(Kg), m_Kf(Kf), m_De(Dg), m_Df(Df)
{
    const auto wa = computeAttachments_();
    const Vector4r qMid = 0.5 * (wa.qA.coeffs() + wa.qB.coeffs());
    const Matrix33r A_mid = Quaternion4r(qMid).normalized().toRotationMatrix();
    m_delta0 = A_mid.transpose() * (wa.xB - wa.xA);

    const Vector4r dQ = wa.qB.coeffs() - wa.qA.coeffs();
    const real_t   factor = (real_t)2.0 / qMid.squaredNorm();
    const real_t   Qmid_w = qMid(3);
    const Vector3r Qmid_q = qMid.head<3>();
    m_kappa0 =  factor * (Qmid_w * dQ.head<3>() - Qmid_q.cross(dQ.head<3>()) - Qmid_q * dQ(3));
}
    

ConstraintResult BeamConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;
    const auto wa = computeAttachments_();

    // Mid-orientation and strains
    const Vector4r qMid = 0.5 * (wa.qA.coeffs() + wa.qB.coeffs());
    const Matrix33r A_mid = Quaternion4r(qMid).normalized().toRotationMatrix();
    const Vector3r gamma = A_mid.transpose() * (wa.xB - wa.xA);   // axial + shear strain
    const Matrix33r gamma_skew = skew_from_vector(gamma);

    const Vector4r dQ = wa.qB.coeffs() - wa.qA.coeffs();
    const real_t   factor = (real_t)2.0 / qMid.squaredNorm();
    const real_t   Qmid_w = qMid(3);
    const Vector3r Qmid_q = qMid.head<3>();
    const Vector3r kappa = factor * (Qmid_w * dQ.head<3>() - Qmid_q.cross(dQ.head<3>()) - Qmid_q * dQ(3));  // bending + twisting strain

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

}} // namespace cardillo::physics
