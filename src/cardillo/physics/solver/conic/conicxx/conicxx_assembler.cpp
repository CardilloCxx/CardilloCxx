#include "conicxx_assembler.hpp"

#include "../conic_scaling.hpp"
#include "../../../assembly/dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

VectorXr ConicxxAssembler::bilateralCompliance(real_t dt, real_t theta) const {
    VectorXr c(m_dyn->numSprings() + m_dyn->numDampers());
    c << m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta), m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta);
    return c;
}

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::P(real_t dt, real_t theta) {
    if (m_compliant_rows) {
        m_P_cache = TripletMatrix::fromDiag(m_dyn->MDiag()).asSparse();
    } else {
        VectorXr diag(m_dyn->numV() + m_dyn->numSprings() + m_dyn->numDampers());
        diag << m_dyn->MDiag(), bilateralCompliance(dt, theta);
        m_P_cache = TripletMatrix::fromDiag(diag).asSparse();
    }
    return m_P_cache;
}

VectorXr& ConicxxAssembler::q(real_t dt, real_t /*theta*/) {
    const int n_pi = m_compliant_rows ? 0 : m_dyn->numSprings() + m_dyn->numDampers();
    m_q_cache.resize(m_dyn->numV() + n_pi);
    m_q_cache << -(m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal())), VectorXr::Zero(n_pi);
    return m_q_cache;
}

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::A(real_t dt, real_t theta) {
    const int n_contact = m_dyn->numContactRows();
    const TripletMatrix contact = m_dyn->W().scaleRows(computeSmu()) * -1.0;
    if (m_compliant_rows) {
        m_A_cache = m_dyn->Wg().vConcat(m_dyn->Wgamma()).vConcat(contact).asSparse();
        return m_A_cache;
    }
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta));
    TripletMatrix A_diag = TripletMatrix::fromDiag(m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta));
    TripletMatrix A = (m_dyn->Wg() | C | TripletMatrix::zero(C.nRows(), A_diag.nCols()))
                          .vConcat(m_dyn->Wgamma() | TripletMatrix::zero(A_diag.nRows(), C.nCols()) | A_diag)
                          .vConcat(contact | TripletMatrix::zero(n_contact, m_dyn->numSprings() + m_dyn->numDampers()));
    m_A_cache = A.asSparse();
    return m_A_cache;
}

VectorXr& ConicxxAssembler::b(real_t dt, real_t theta) {
    const auto& Wcontact = m_dyn->W().asSparse();
    VectorXr Smu = computeSmu();
    VectorXr restitution = m_dyn->restitutionVec();
    VectorXr biasImpulse = m_dyn->contactVVec();

    VectorXr gamma_old = Wcontact * m_dyn->vVec() + m_dyn->contactVVec();

    if ((restitution.array() != 0.0).any())
        biasImpulse += restitution.cwiseProduct(gamma_old);

    VectorXr b_top = -m_dyn->springBias(dt, theta);
    VectorXr b_mid = -m_dyn->damperBias(dt, theta);
    VectorXr b_contact = scaledContactOffset(*m_dyn, Smu, biasImpulse, gamma_old);

    m_b_cache.resize(b_top.size() + b_mid.size() + b_contact.size());
    m_b_cache << b_top, b_mid, b_contact;

    return m_b_cache;
}

const VectorXr& ConicxxAssembler::compliance(real_t dt, real_t theta) {
    if (m_compliant_rows) {
        m_compliance_cache = bilateralCompliance(dt, theta);
    } else {
        m_compliance_cache.setZero(m_dyn->numSprings() + m_dyn->numDampers());
    }
    return m_compliance_cache;
}

const conicxx::ConeSpec& ConicxxAssembler::coneSpec() {
    m_cone_spec_cache = conicxx::ConeSpec{};
    m_cone_spec_cache.zero_dim = m_dyn->numSprings() + m_dyn->numDampers();
    m_cone_spec_cache.nonneg_dim = m_dyn->numFrictionlessContacts();

    const int nsoc = m_dyn->numFrictionalContacts();
    m_cone_spec_cache.soc_dims.assign(static_cast<std::size_t>(nsoc), 3);

    return m_cone_spec_cache;
}

VectorXr ConicxxAssembler::computeSmu() const { return frictionConeScaling(*m_dyn); }

}  // namespace cardillo::physics::assembly
