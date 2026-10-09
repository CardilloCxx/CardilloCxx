#include "clarabel_assembler.hpp"

#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

const SparseMatrix<Eigen::ColMajor>& ClarabelAssembler::P(real_t dt, real_t theta) {
    auto M = m_dyn->MDiag();
    auto C = m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta);
    auto A = m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta);

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    TripletMatrix P = TripletMatrix::fromDiag(diag);
    m_P_cache = P.asSparse();
    return m_P_cache;
}

VectorXr& ClarabelAssembler::q(real_t dt, real_t theta) {
    (void)theta;

    auto q_top = -(m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal()));
    auto q_bot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    m_q_cache.resize(q_top.size() + q_bot.size());
    m_q_cache << q_top, q_bot;

    return m_q_cache;
}

const SparseMatrix<Eigen::ColMajor>& ClarabelAssembler::A(real_t dt, real_t theta) {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta));
    TripletMatrix A_diag = TripletMatrix::fromDiag(m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta));

    const int n_contact = m_dyn->numContactRows();
    VectorXr Smu = computeSmu();

    TripletMatrix A = (m_dyn->Wg() | (C * -1.0) | TripletMatrix::zero(C.nRows(), A_diag.nCols()))
                          .vConcat(m_dyn->Wgamma() | TripletMatrix::zero(A_diag.nRows(), C.nCols()) | (A_diag * -1.0))
                          .vConcat(m_dyn->W().scaleRows(Smu) * -1.0 | TripletMatrix::zero(n_contact, m_dyn->numSprings() + m_dyn->numDampers()));

    m_A_cache = A.asSparse();
    return m_A_cache;
}

VectorXr& ClarabelAssembler::b(real_t dt, real_t theta) {
    const auto& Wcontact = m_dyn->W().asSparse();
    VectorXr Smu = computeSmu();
    VectorXr restitution = m_dyn->restitutionVec();
    VectorXr biasImpulse = m_dyn->contactVVec();

    VectorXr gamma_old = Wcontact * m_dyn->vVec() + m_dyn->contactVVec();

    if ((restitution.array() != 0.0).any())
        biasImpulse += restitution.cwiseProduct(gamma_old);

    VectorXr b_top = -m_dyn->springBias(dt, theta);
    VectorXr b_mid = -m_dyn->damperBias(dt, theta);
    VectorXr b_contact = Smu.cwiseProduct(biasImpulse);

    for (int i = m_dyn->numFrictionlessContacts(); i < m_dyn->numContactRows(); i += 3) {
        const real_t y_1 = gamma_old[i + 1];
        const real_t y_2 = gamma_old[i + 2];
        b_contact[i] += std::sqrt(y_1 * y_1 + y_2 * y_2);
    }

    m_b_cache.resize(b_top.size() + b_mid.size() + b_contact.size());
    m_b_cache << b_top, b_mid, b_contact;

    return m_b_cache;
}

const std::vector<clarabel::SupportedConeT<double>>& ClarabelAssembler::cones() {
    m_cones_cache.clear();

    const uintptr_t p = static_cast<uintptr_t>(m_dyn->numSprings() + m_dyn->numDampers());
    const uintptr_t l = static_cast<uintptr_t>(m_dyn->numFrictionlessContacts());
    const uintptr_t nsoc = static_cast<uintptr_t>(m_dyn->numFrictionalContacts());

    if (p > 0) {
        m_cones_cache.push_back(clarabel::ZeroConeT<double>(p));
    }
    if (l > 0) {
        m_cones_cache.push_back(clarabel::NonnegativeConeT<double>(l));
    }
    for (uintptr_t i = 0; i < nsoc; ++i) {
        (void)i;
        m_cones_cache.push_back(clarabel::SecondOrderConeT<double>((uintptr_t)3));
    }

    return m_cones_cache;
}

VectorXr ClarabelAssembler::computeSmu() const {
    VectorXr Smu = VectorXr::Ones(m_dyn->numContactRows());
    for (const auto& c : m_dyn->contacts()) {
        if (c.friction_mu <= 0) continue;
        const int idx = c.impulse_base_index;
        if (idx < 0 || idx >= Smu.size()) continue;
        Smu[idx] = 1.0 / c.friction_mu;
    }
    return Smu;
}

}  // namespace cardillo::physics::assembly
