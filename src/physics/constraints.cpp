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

// ===================== TranslationalConstraint =====================
ConstraintResult TranslationalConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;
    const auto wa = computeAttachments_();

//         out.WgA = MatrixXXr::Zero(6, 3);
//         out.WgA.block<3,3>(0,0) = -wa.RA;
//         out.WgA.block<3,3>(3,0) = -skew_from_vector(m_rA_local);
// 
//         out.WgB = MatrixXXr::Zero(6, 3);
//         out.WgB.block<3,3>(0,0) = wa.RA;
//         out.WgB.block<3,3>(3,0) = (skew_from_vector(m_rB_local) * wa.RB.transpose() * wa.RA);

    out.WgA = MatrixXXr::Zero(6, 3);
    out.WgA.block<3,3>(0,0) = -wa.RA;
    out.WgA.block<3,3>(3,0) = -skew_from_vector(wa.RA.transpose() * (wa.xB - wa.xA) + m_rA_local);

    out.WgB = MatrixXXr::Zero(6, 3);
    out.WgB.block<3,3>(0,0) = wa.RA;
    out.WgB.block<3,3>(3,0) = skew_from_vector(m_rB_local) * wa.RB.transpose() * wa.RA;

    // (-wa.RA.transpose() * wa.RB * skew_from_vector(m_rB_local)).transpose();
    // out.WgA.block<3,3>(3,0) = (skew_from_vector(wa.RA.transpose() * (wa.xB - wa.xA)) + wa.RA.transpose() * skew_from_vector(ra_world) * wa.RA).transpose();
    // out.WgB.block<3,3>(3,0) = (wa.RA.transpose() * skew_from_vector(rb_world) * wa.RB).transpose();

    out.WgammaA = out.WgA;
    out.WgammaB = out.WgB;

    out.Crows = VectorXr::Zero(3);
    out.Arows = VectorXr::Zero(3);
    fillCompliance3(out.Crows, 0, m_K);
    fillCompliance3(out.Arows, 0, m_D);

    return out;
}

// ===================== RigidConstraint =====================
ConstraintResult RigidConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;
    const auto wa = computeAttachments_();

    // Prepare 6 scalar rows (3 translational + 3 rotational) with 6 columns per velocity block
    out.WgA = MatrixXXr::Zero(6, 6);
    out.WgB = MatrixXXr::Zero(6, 6);

    // Lock translations in A's frame (columns 0..2)
    out.WgA.block<3,3>(0,0) = -wa.RA;
    out.WgA.block<3,3>(3,0) = -skew_from_vector(wa.RA.transpose() * (wa.xB - wa.xA) + m_rA_local);
    out.WgB.block<3,3>(0,0) = wa.RA;
    out.WgB.block<3,3>(3,0) = skew_from_vector(m_rB_local) * wa.RB.transpose() * wa.RA;

    // Lock rotations in A's frame (columns 3..5)
    out.WgA.block<3,3>(3,3) = -Matrix33r::Identity();
    out.WgB.block<3,3>(3,3) = Matrix33r::Identity();

    out.WgammaA = out.WgA;
    out.WgammaB = out.WgB;

    out.Crows = VectorXr::Constant(6, std::numeric_limits<real_t>::infinity());
    out.Arows = VectorXr::Zero(6);

    return out;
}

// ===================== HingeConstraint =====================
HingeConstraint::HingeConstraint(entt::registry& reg,
                                 entt::entity a,
                                 entt::entity b,
                                 const Vector3r& rA_local,
                                 const Vector3r& rB_local,
                                 const Vector3r& axis_localA,
                                 const real_t& K_axis,
                                 const real_t& D_axis,
                                 const Vector2r& Kf_A,
                                 const Vector2r& Df_A,
                                 const Vector3r& Ke_A,
                                 const Vector3r& De_A)
    : ConstraintPattern(reg, a, b, rA_local, rB_local),
      m_hingeFrame(Matrix33r::Identity()),
      m_Ke(Ke_A), m_De(De_A),
      m_Kf(Vector3r(Kf_A.x(), Kf_A.y(), K_axis)),
      m_Df(Vector3r(Df_A.x(), Df_A.y(), D_axis))
{
    // Compute hinge frame in A's local space
    Vector3r axisA = axis_localA.normalized();
    Vector3r up = Vector3r(0, 0, 1);
    if (std::abs(axisA.dot(up)) > (real_t)0.99) { up = Vector3r(0, 1, 0); }
    Vector3r right = up.cross(axisA).normalized();
    up = -right.cross(axisA);
    m_hingeFrame.col(0) = right;
    m_hingeFrame.col(1) = up;
    m_hingeFrame.col(2) = axisA;
}
ConstraintResult HingeConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;
    const auto wa = computeAttachments_();

    // Prepare 6 scalar rows (3 translational + 3 rotational) with 6 columns per velocity block
    out.WgA = MatrixXXr::Zero(6, 6);
    out.WgB = MatrixXXr::Zero(6, 6);

    // Lock translations in A's frame (columns 0..2)
    out.WgA.block<3,3>(0,0) = -wa.RA;
    out.WgA.block<3,3>(3,0) = -skew_from_vector(wa.RA.transpose() * (wa.xB - wa.xA) + m_rA_local);

    out.WgB.block<3,3>(0,0) = wa.RA;
    out.WgB.block<3,3>(3,0) = skew_from_vector(m_rB_local) * wa.RB.transpose() * wa.RA;

    // Lock rotation relative in A's local frame
    out.WgA.block<3,3>(3,3) = -m_hingeFrame;
    out.WgB.block<3,3>(3,3) =  m_hingeFrame * wa.RA.transpose() *  wa.RB;

    out.WgammaA = out.WgA;
    out.WgammaB = out.WgB;

    out.Crows = VectorXr::Zero(6);
    out.Arows = VectorXr::Zero(6);
    fillCompliance3(out.Crows, 0, m_Ke);
    fillCompliance3(out.Crows, 3, m_Kf);
    fillCompliance3(out.Arows, 0, m_De);
    fillCompliance3(out.Arows, 3, m_Df);

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
