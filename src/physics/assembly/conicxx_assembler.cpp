#include "conicxx_assembler.hpp"

#include "conic_scaling.hpp"
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::P(real_t dt, real_t theta) {
    auto M = m_dyn->MDiag();
    auto C = m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta);
    auto A = m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta);

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    TripletMatrix P = TripletMatrix::fromDiag(diag);
    m_P_cache = P.asSparse();
    return m_P_cache;
}

VectorXr& ConicxxAssembler::q(real_t dt, real_t /*theta*/) {

    auto q_top = -(m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal()));
    auto q_bot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    m_q_cache.resize(q_top.size() + q_bot.size());
    m_q_cache << q_top, q_bot;

    return m_q_cache;
}

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::A(real_t dt, real_t theta) {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta));
    TripletMatrix A_diag = TripletMatrix::fromDiag(m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta));

    const int n_contact = m_dyn->numContactRows();
    VectorXr Smu = computeSmu();

    TripletMatrix A = (m_dyn->Wg() | C | TripletMatrix::zero(C.nRows(), A_diag.nCols()))
                          .vConcat(m_dyn->Wgamma() | TripletMatrix::zero(A_diag.nRows(), C.nCols()) | A_diag)
                          .vConcat(m_dyn->W().scaleRows(Smu) * -1.0 | TripletMatrix::zero(n_contact, m_dyn->numSprings() + m_dyn->numDampers()));

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
