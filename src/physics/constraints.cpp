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
                                                             const JointFrame& frame,
                                                             const Vector3r& K_trans,
                                                             const Vector3r& D_trans,
                                                             const Vector3r& K_rot,
                                                             const Vector3r& D_rot)
    : ConstraintPattern(reg, a, b,Vector3r::Zero(), Vector3r::Zero())
    , m_joint(JointProperties(reg, a, b, frame))
    , m_K_trans(K_trans)
    , m_D_trans(D_trans)
    , m_K_rot(K_rot)
    , m_D_rot(D_rot)
{
    m_rA_local = m_joint.K1_r_S1J;
    m_rB_local = m_joint.K2_r_S2J;
}

// Helper: build full 6x6 Jacobians for a joint between A and B using
// the precomputed joint-frame geometry. This encodes a
// fully locked 6-DOF joint; specialised constraints will mask DOFs via C/A.
void TranslationRotationConstraint::buildJointJacobian(const ConstraintPattern::WorldAttachments& wa,
                                                       MatrixXXr& WgA,
                                                       MatrixXXr& WgB) const {

    WgA = MatrixXXr::Zero(6, 6);
    WgB = MatrixXXr::Zero(6, 6);

    const auto& A_IK1 = wa.RA;
    const auto& A_IK2 = wa.RB;
    const auto& A_K1J = m_joint.A_K1J;
    const Matrix33r A_IJ = A_IK1 * A_K1J;
    const Matrix33r A_K2J =  A_IK2.transpose() * A_IJ;
    const Matrix33r skew_g = skew_from_vector(m_joint.compute_g(wa.pA, wa.pB, wa.RA, wa.RB));
    
    // translations
    WgA.block<3,3>(0,0) = -A_IJ; 
    WgA.block<3,3>(3,0) = -m_joint.K1_r_S1J_skew * A_K1J - A_K1J * skew_g;
    WgB.block<3,3>(0,0) = A_IJ;
    WgB.block<3,3>(3,0) = m_joint.K2_r_S2J_skew * A_K2J;

    // orientations
    WgA.block<3,3>(3,3) = A_K1J;
    WgB.block<3,3>(3,3) = -A_K2J;
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
BeamConstraint::BeamConstraint( entt::registry& reg,
                                entt::entity a,
                                entt::entity b,
                                const cardillo::PhysicsSystem::BeamSpringParams& springs,
                                const cardillo::PhysicsSystem::BeamCrossSection& section)
        : ConstraintPattern(reg, a, b, Vector3r::Zero(), Vector3r::Zero()),
            m_springs(springs), m_section(section)
{
    m_crackStrainMax = springs.crackStrainMax;
    const auto wa = computeAttachments_();
    const Vector4r qMid = 0.5 * (wa.qA.coeffs() + wa.qB.coeffs());
    const Matrix33r A_mid = Quaternion4r(qMid).normalized().toRotationMatrix();
    m_delta0 = A_mid.transpose() * (wa.xB - wa.xA);
    l_0 = m_delta0.norm();

    const Vector4r dQ = wa.qB.coeffs() - wa.qA.coeffs();
    const real_t   factor = (real_t)2.0 / qMid.squaredNorm();
    const real_t   Qmid_w = qMid(3);
    const Vector3r Qmid_q = qMid.head<3>();
    m_kappa0 =  factor * (Qmid_w * dQ.head<3>() - Qmid_q.cross(dQ.head<3>()) - Qmid_q * dQ(3));
}

real_t CalculateMaxPrincipalTensileStress(
    const Vector3r& gamma, 
    const Vector3r& kappa,
    const Vector3r& Ke_eff,
    const Vector3r& Kf_eff,
    const PhysicsSystem::BeamCrossSection& sec
) {
    // Axial Force (N)
    real_t N  = Ke_eff.x() * gamma.x();
    
    // Bending Moments (My, Mz)
    real_t My = Kf_eff.y() * kappa.y();
    real_t Mz = Kf_eff.z() * kappa.z();
    
    // Shear Forces (Vy, Vz)
    real_t Vy = Ke_eff.y() * gamma.y();
    real_t Vz = Ke_eff.z() * gamma.z();
    
    const real_t A = sec.area();
    const real_t W = sec.sectionModulus(); // Assume sec.sectionModulus() returns I/R 
    
    // Safety check: Avoid division by zero
    if (A <= (real_t)0.0 || W <= (real_t)0.0) return (real_t)0.0;
    
    // Axial Stress Component (N/A)
    real_t sigma_axial = N / A; 
    
    // Combined Bending Moment Magnitude: Mb = sqrt(My^2 + Mz^2)
    real_t Mb_mag = std::sqrt(My * My + Mz * Mz);
    
    // Maximum Tensile Axial/Bending Stress (Sigma_xx) at the perimeter:
    real_t sigma_xx_max = sigma_axial + Mb_mag / W; 
    
    // Shear Stress (tau_xy, tau_xz) at this point:
    // For a circular section, shear stress is ZERO at the perimeter where bending is max.
    // The Rankine criterion simplifies significantly here:
    real_t tau_xy_at_max = 0.0;
    real_t tau_xz_at_max = 0.0;
    
    // Sigma_1 = (sigma_xx / 2) + sqrt((sigma_xx / 2)^2 + tau_xy^2 + tau_xz^2)
    // Since tau_xy = tau_xz = 0, this simplifies to Sigma_1 = Sigma_xx_max (if tensile)
    real_t sigma_1 = (sigma_xx_max / 2.0) + std::sqrt(
        std::pow(sigma_xx_max / 2.0, 2) + 
        std::pow(tau_xy_at_max, 2) + 
        std::pow(tau_xz_at_max, 2)
    );
    
    // We only care about TENSION for the Rankine failure criterion
    return std::max((real_t)0.0, sigma_1);
}
    

ConstraintResult BeamConstraint::getConstraint() const {
    ConstraintResult out; out.a = m_a; out.b = m_b;
    const auto wa = computeAttachments_();

    // Mid-orientation and strains (robust): enforce quaternion hemisphere continuity and use slerp
    Quaternion4r qA = wa.qA;
    Quaternion4r qB = wa.qB;
    const Quaternion4r qMidQ(0.5 * (qA.coeffs() + qB.coeffs()));   // 0.5*(qA + qB) 
    const Matrix33r A_mid = qMidQ.normalized().toRotationMatrix();
    const Vector3r gamma = A_mid.transpose() * (wa.xB - wa.xA);   // axial + shear strain
    const Matrix33r gamma_skew = skew_from_vector(gamma);

    const Vector4r dQ = qB.coeffs() - qA.coeffs();
    const real_t   factor = (real_t)2.0 / qMidQ.coeffs().squaredNorm();
    const real_t   Qmid_w = qMidQ.coeffs()(3);
    const Vector3r Qmid_q = qMidQ.coeffs().head<3>();
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

    // Effective stiffnesses (considering scaling factors)
    Vector3r Ke_eff = m_springs.Ke(l_0, m_section).cwiseProduct(m_springs.scaleKe);
    Vector3r Kf_eff = m_springs.Kf(l_0, m_section).cwiseProduct(m_springs.scaleKf);
    Vector3r De_eff = m_springs.De;
    Vector3r Df_eff = m_springs.Df;

    if (!std::isfinite( m_springs.tensileStrength)) { 
        fillCompliance3(out.Crows, 0, Ke_eff);
        fillCompliance3(out.Crows, 3, Kf_eff);
        fillCompliance3(out.Arows, 0, De_eff);
        fillCompliance3(out.Arows, 3, Df_eff);
        return out;
    }

    // 1. Calculate Stresses (N, My, Mz, Vy, Vz) -> then stress tensor and sigma_1 at critical point(s)
    real_t sigma_1 = CalculateMaxPrincipalTensileStress(
        gamma, kappa,
        Ke_eff,
        Kf_eff,
        m_section
    );

    const Vector3r strain_e = gamma - m_delta0; 
    real_t e_ck_current = std::max(std::abs(strain_e.y() / l_0), std::abs(strain_e.z() / l_0));

    // 2. Rankine Criterion for Crack Initiation
    if (!m_crackReported && sigma_1 > m_springs.tensileStrength) {
        m_crackReported = true;
        m_crackStrainPeak = e_ck_current; // Initialize the peak strain when crack first forms
    }

    // 3. Shear Retention Driven by Peak Crack Strain (post-initiation)
    real_t retention = (real_t)1.0;
    if (m_crackReported) {
        // Only update the peak strain AFTER the crack has been initiated by Rankine
        if (e_ck_current > m_crackStrainPeak) m_crackStrainPeak = e_ck_current;
        
        // Apply the softening law (Retention = 1 - e_ck_peak / e_ck_max)
        if (m_crackStrainPeak > 0) {
            retention = std::clamp((real_t)1.0 - (m_crackStrainPeak / m_springs.crackStrainMax), (real_t)0.0, (real_t)1.0);
        }
    }

    // 4. Apply Retention and Check for Permanent Break (The rest of your code)
    Ke_eff.y() *= retention;
    Ke_eff.z() *= retention;
    // Break permanently once retention reaches 0
    if (retention <= (real_t)0.0) m_broken = true;

    if (m_broken) {
        if (m_reg && m_reg->valid(m_a) && m_reg->all_of<cardillo::PhysicsSystem::C_BeamElement>(m_a)) {
            auto &beam_a = m_reg->get<cardillo::PhysicsSystem::C_BeamElement>(m_a);
            beam_a.next = std::nullopt;
            if (!beam_a.next.has_value() && !beam_a.prev.has_value()) {
                m_reg->remove<cardillo::PhysicsSystem::C_BeamElement>(m_a);
            }
        }
        if (m_reg && m_reg->valid(m_b) && m_reg->all_of<cardillo::PhysicsSystem::C_BeamElement>(m_b)) {
            auto &beam_b = m_reg->get<cardillo::PhysicsSystem::C_BeamElement>(m_b);
            beam_b.prev = std::nullopt;
            if (!beam_b.next.has_value() && !beam_b.prev.has_value()) {
                m_reg->remove<cardillo::PhysicsSystem::C_BeamElement>(m_b);
            }
        }
        Ke_eff = Vector3r::Zero();
        Kf_eff = Vector3r::Zero();
        De_eff = Vector3r::Zero();
        Df_eff = Vector3r::Zero();
    }

    fillCompliance3(out.Crows, 0, Ke_eff);
    fillCompliance3(out.Crows, 3, Kf_eff);
    fillCompliance3(out.Arows, 0, De_eff);
    fillCompliance3(out.Arows, 3, Df_eff);

    return out;
}

}} // namespace cardillo::physics
