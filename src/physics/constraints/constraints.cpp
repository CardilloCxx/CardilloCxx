#include "constraints.hpp"

namespace cardillo {
namespace physics {

// ===================== ConstraintPattern base =====================
ConstraintPattern::ConstraintPattern(entt::registry& reg, entt::entity a, entt::entity b, index_t numRows) : m_reg(&reg), m_a(a), m_b(b), m_numRows(numRows) {}

bool ConstraintPattern::getAttachPointsWorld(Vector3r& xA, Vector3r& xB) const {
    auto wa = computeAttachments_();
    xA = wa.xA;
    xB = wa.xB;
    return true;
}

ConstraintPattern::WorldAttachments ConstraintPattern::computeAttachments_() const {
    WorldAttachments wa;
    if (m_reg) {
        const auto stateA = cardillo::RigidBody::getState(*m_reg, m_a);
        const auto stateB = cardillo::RigidBody::getState(*m_reg, m_b);

        wa.qA = stateA.orientation;
        wa.qB = MathHelper::alignQuaternionTo(stateB.orientation, wa.qA);
        wa.RA = stateA.rotation;
        wa.RB = stateB.rotation;
        wa.rA_world = wa.RA * m_rA_local;
        wa.rB_world = wa.RB * m_rB_local;
        wa.pA = stateA.position;
        wa.pB = stateB.position;
        wa.xA = wa.pA + wa.rA_world;
        wa.xB = wa.pB + wa.rB_world;

        if (wa.qA.coeffs().dot(wa.qB.coeffs()) < (real_t)0.0) {
            wa.qB.coeffs() = -wa.qB.coeffs();
        }
    }
    return wa;
}

real_t ConstraintPattern::stiffnessToCompliance(real_t k) {
    const real_t eps = (real_t)1e-12;
    if (std::isfinite(k)) return (k > eps) ? (real_t)1 / k : std::numeric_limits<real_t>::infinity();
    return (real_t)0;  // +inf stiffness -> zero compliance
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
    ConstraintResult out;
    out.a = m_a;
    out.b = m_b;

    const auto wa = computeAttachments_();
    const Vector3r r = wa.xB - wa.xA;
    const real_t len = r.norm();
    out.positionError = VectorXr::Constant(1, len - L0);

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
    out.WgB.block(0, 0, 3, 1) = n;
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

TranslationRotationConstraint::TranslationRotationConstraint(entt::registry& reg, entt::entity a, entt::entity b, const JointFrame& frame, const Vector3r& K_trans, const Vector3r& D_trans,
                                                             const Vector3r& K_rot, const Vector3r& D_rot)
    : ConstraintPattern(reg, a, b, 6), m_joint(JointProperties(reg, a, b, frame)), m_K_trans(K_trans), m_D_trans(D_trans), m_K_rot(K_rot), m_D_rot(D_rot) {
    m_rA_local = m_joint.K1_r_S1J;
    m_rB_local = m_joint.K2_r_S2J;

    const auto res = getConstraint();
    m_g0 = res.positionError;
}

// Helper: build full 6x6 Jacobians for a joint between A and B using
// the precomputed joint-frame geometry. This encodes a
// fully locked 6-DOF joint; specialised constraints will mask DOFs via C/A.
void TranslationRotationConstraint::buildJointJacobian(const ConstraintPattern::WorldAttachments& wa, Vector3r g, MatrixXXr& WgA, MatrixXXr& WgB) const {
    WgA = MatrixXXr::Zero(6, 6);
    WgB = MatrixXXr::Zero(6, 6);

    const auto& A_IK1 = wa.RA;
    const auto& A_IK2 = wa.RB;
    const auto& A_K1J = m_joint.A_K1J;
    const Matrix33r A_IJ = A_IK1 * A_K1J;
    const Matrix33r A_K2J = A_IK2.transpose() * A_IJ;
    const Matrix33r skew_g = skew_from_vector(g);

    // translations
    WgA.block<3, 3>(0, 0) = -A_IJ;
    WgA.block<3, 3>(3, 0) = -m_joint.K1_r_S1J_skew * A_K1J - A_K1J * skew_g;
    WgB.block<3, 3>(0, 0) = A_IJ;
    WgB.block<3, 3>(3, 0) = m_joint.K2_r_S2J_skew * A_K2J;

    // orientations
    WgA.block<3, 3>(3, 3) = A_K1J;
    WgB.block<3, 3>(3, 3) = -A_K2J;
}

ConstraintResult TranslationRotationConstraint::getConstraint() const {
    ConstraintResult out;
    out.a = m_a;
    out.b = m_b;

    // Use standard attachment computation; joint is defined by m_rA_local/m_rB_local
    const auto wa = computeAttachments_();

    // Build full 6x6 Jacobians for a rigid joint at the current attachments.
    const Vector3r g = m_joint.compute_g(wa.pA, wa.pB, wa.RA, wa.RB);
    buildJointJacobian(wa, g, out.WgA, out.WgB);

    out.positionError = getPositionError(g, out);

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

VectorXr TranslationRotationConstraint::getPositionError(Vector3r g, ConstraintResult res) const {
    VectorXr posErr(6);
    posErr.head<3>() = g;
    const auto R = -(res.WgA.block<3, 3>(3, 3)).transpose() * res.WgB.block<3, 3>(3, 3);
    posErr.tail<3>() = Vector3r(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1)) * (real_t)0.5;
    return posErr - m_g0;
}

// ===================== BeamConstraint =====================
BeamConstraint::BeamConstraint(entt::registry& reg, entt::entity a, entt::entity b, const cardillo::physics::BeamSpringParams& springs, const cardillo::physics::BeamCrossSection& section)
    : ConstraintPattern(reg, a, b, 6), m_springs(springs), m_section(section) {
    const auto res = getConstraint();
    m_gamma0 = res.positionError.head<3>();
    m_kappa0 = res.positionError.tail<3>();
}

ConstraintResult BeamConstraint::getConstraint() const {
    ConstraintResult out;
    out.a = m_a;
    out.b = m_b;
    const auto wa = computeAttachments_();

    // Mid-orientation and strains (robust): enforce quaternion hemisphere continuity and use slerp
    const Quaternion4r qMidQ(0.5 * (wa.qA.coeffs() + wa.qB.coeffs()));  // 0.5*(qA + qB)
    const Matrix33r A_mid = qMidQ.normalized().toRotationMatrix();

    const Vector3r gamma = A_mid.transpose() * (wa.xB - wa.xA);  // axial + shear strain

    const Vector4r dQ = wa.qB.coeffs() - wa.qA.coeffs();
    const real_t factor = (real_t)2.0 / qMidQ.coeffs().squaredNorm();
    const real_t Qmid_w = qMidQ.coeffs()(3);
    const Vector3r Qmid_q = qMidQ.coeffs().head<3>();
    const Vector3r kappa = factor * (Qmid_w * dQ.head<3>() - Qmid_q.cross(dQ.head<3>()) - Qmid_q * dQ(3));  // bending + twisting strain

    out.positionError = VectorXr::Zero(6);
    out.positionError.head<3>() = gamma - m_gamma0;
    out.positionError.tail<3>() = kappa - m_kappa0;

    const Matrix33r gamma_skew = skew_from_vector(gamma);
    const Matrix33r kappa_skew = skew_from_vector(kappa);

    out.WgA = MatrixXXr::Zero(6, 6);
    out.WgA.block<3, 3>(0, 0) = -A_mid;
    out.WgA.block<3, 3>(3, 0) = (real_t)-0.5 * gamma_skew;
    out.WgA.block<3, 3>(3, 3) = -Matrix33r::Identity() - (real_t)0.5 * kappa_skew;  // l_0 jeweils rausgekürzt

    out.WgB = MatrixXXr::Zero(6, 6);
    out.WgB.block<3, 3>(0, 0) = A_mid;
    out.WgB.block<3, 3>(3, 0) = (real_t)-0.5 * gamma_skew;
    out.WgB.block<3, 3>(3, 3) = Matrix33r::Identity() - (real_t)0.5 * kappa_skew;  // l_0 jeweils rausgekürzt

    out.WgammaA = out.WgA;
    out.WgammaB = out.WgB;

    // Row-wise compliances (spring and damper)
    out.Crows = VectorXr::Zero(6);
    out.Arows = VectorXr::Zero(6);

    // Effective stiffnesses (considering scaling factors)
    Vector3r Ke_eff = m_springs.Ke(l_0, m_section).cwiseProduct(m_springs.scaleKe);
    Vector3r Kf_eff = m_springs.Kf(l_0, m_section).cwiseProduct(m_springs.scaleKf);
    Vector3r De_eff = Ke_eff * m_springs.dampingFactor;
    Vector3r Df_eff = Kf_eff * m_springs.dampingFactor;

    fillCompliance3(out.Crows, 0, Ke_eff);
    fillCompliance3(out.Crows, 3, Kf_eff);
    fillCompliance3(out.Arows, 0, De_eff);
    fillCompliance3(out.Arows, 3, Df_eff);

    return out;
}
}  // namespace physics
}  // namespace cardillo
