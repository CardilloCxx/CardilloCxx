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

static Vector3r so3_log_from_R(const Matrix33r &R) {
    // Use trace to compute angle
    real_t cos_th = (R.trace() - 1.0) * 0.5;
    cos_th = std::min((real_t)1.0, std::max((real_t)-1.0, cos_th));
    real_t th = std::acos(cos_th);

    if (th < 1e-8) {
        // small: phi ≈ 0.5 * vee(R - R^T)
        Vector3r v;
        v.x() = (R(2,1) - R(1,2)) * 0.5;
        v.y() = (R(0,2) - R(2,0)) * 0.5;
        v.z() = (R(1,0) - R(0,1)) * 0.5;
        return v;
    } else {
        real_t sin_th = std::sin(th);
        Vector3r axis;
        axis.x() = (R(2,1) - R(1,2)) / (2.0 * sin_th);
        axis.y() = (R(0,2) - R(2,0)) / (2.0 * sin_th);
        axis.z() = (R(1,0) - R(0,1)) / (2.0 * sin_th);
        return axis * th;
    }
}

// left-Jacobian inverse J^{-1}(phi) with series for small phi
static Matrix33r leftJacobianInverse(const Vector3r &phi) {
    const real_t th = phi.norm();
    const Matrix33r I = Matrix33r::Identity();
    const Matrix33r Phi = skew_from_vector(phi);

    if (th < 1e-8) {
        // Series: J^{-1} ≈ I + 0.5 Phi + (1/12) Phi^2
        return I + (real_t)0.5 * Phi + (real_t)(1.0/12.0) * (Phi * Phi);
    } else {
        const real_t half_th = 0.5 * th;
        const real_t cot_half = 1.0 / std::tan(half_th);
        // exact formula: J^{-1} = I + 0.5 Phi + (1/th^2) * (1 - (th * cot(th/2))/2) * Phi^2
        const real_t th2 = th * th;
        const real_t factor = (1.0 / th2) * (1.0 - (th * cot_half) * 0.5);
        return I + (real_t)0.5 * Phi + factor * (Phi * Phi);
    }
}

// ===================== RigidConstraint =====================
ConstraintResult RigidConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;

    // TODO: Move this in computeAttachments_ or a similar function
    const Vector3r& r_OS1 = m_reg->get<cardillo::PhysicsSystem::C_Position3>(m_a).value;
    const Vector3r& r_OS2 = m_reg->get<cardillo::PhysicsSystem::C_Position3>(m_b).value;
    const Quaternion4r& Q1 = m_reg->get<cardillo::PhysicsSystem::C_Orientation>(m_a).value;
    const Quaternion4r& Q2 = m_reg->get<cardillo::PhysicsSystem::C_Orientation>(m_b).value;
    const Matrix33r A_IK1 = Q1.toRotationMatrix();
    const Matrix33r A_IK2 = Q2.toRotationMatrix();

    out.WgA = MatrixXXr::Zero(6, 6);
    out.WgB = MatrixXXr::Zero(6, 6);
    
    // translations
    out.WgA.block<3,3>(0,0) = -Matrix33r::Identity();
    out.WgA.block<3,3>(3,0) = -m_K1_r_S1J_skew * A_IK1.transpose();

    out.WgB.block<3,3>(0,0) = Matrix33r::Identity();
    out.WgB.block<3,3>(3,0) = m_K2_r_S2J_skew * A_IK2.transpose();

    // orientations
    // 1. x x y = z
    out.WgA.block<3,1>(3,3) = m_A_K1J.col(2);
    out.WgB.block<3,1>(3,3) = -m_A_K2J.col(2);

    // 2. y @ z = x
    out.WgA.block<3,1>(3,4) = m_A_K1J.col(0);
    out.WgB.block<3,1>(3,4) = -m_A_K2J.col(0);

    // 3. z x x = y
    out.WgA.block<3,1>(3,5) = m_A_K1J.col(1);
    out.WgB.block<3,1>(3,5) = -m_A_K2J.col(1);

    // possibly set compilance
    // TODO: Can be done once since Crows and Arows are constant!
    // TODO: Remove parts of W and C for zero stiffnesses!
    out.Crows = VectorXr::Zero(6);
    if (abs(m_k_axis(0)) > 0) out.Crows(3) = 1.0 / m_k_axis(0);
    if (abs(m_k_axis(1)) > 0) out.Crows(4) = 1.0 / m_k_axis(1);
    if (abs(m_k_axis(2)) > 0) out.Crows(5) = 1.0 / m_k_axis(2);

    // possibly set attenuation
    const Array3b nonzero_d = m_d_axis.array().abs() > 0;
    const index_t nD = nonzero_d.count();

    out.Arows.resize(nD);
    out.Arows = VectorXr::Constant(nD, std::numeric_limits<real_t>::max());
    // TODO: This logic is still missing
    // out.WgammaA.resize(6, nD);
    // out.WgammaB.resize(6, nD);
    // index_t offset = 0;
    // for (index_t i=0; i < 3; ++i) {
    //     if (nonzero_d(i)) {
    //         out.Arows(offset) = 1.0 / m_d_axis(i);
    //         out.WgammaA.col(offset) = out.WgA.col(3 + i);
    //         out.WgammaB.col(offset) = out.WgA.col(3 + i);
    //         offset++;
    //     }
    // }

    return out;
}

// ===================== HingeConstraint =====================
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

    // // Lock rotation relative in A's local frame
    out.WgA.block<3,3>(3,3) = -m_hingeFrame;
    out.WgB.block<3,3>(3,3) =  m_hingeFrame * wa.RA.transpose() *  wa.RB;

//     Matrix33r hingeFrame = wa.RA * m_hingeFrame; // hinge frame in world
//     Vector3r hingeAxis = hingeFrame.col(2);      // z-axis of hinge frame
//     Matrix33r P = Matrix33r::Identity() - hingeAxis * hingeAxis.transpose();
// 
//     const Matrix33r Rrel = wa.RA.transpose() * wa.RB;
//     const Vector3r phi = so3_log_from_R(Rrel);
//     const Matrix33r Jinv = leftJacobianInverse(phi);
//     
//     // Lock rotations in A's frame (columns 3..5)
//     out.WgA.block<3,3>(3,3) = -P * Jinv.transpose();
//     out.WgB.block<3,3>(3,3) =  P * (Jinv * Rrel).transpose();

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
