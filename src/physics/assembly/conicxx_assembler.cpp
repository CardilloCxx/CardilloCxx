#include "conicxx_assembler.hpp"

#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::P(real_t dt, real_t theta) {
    auto M = m_dyn->MDiag();
    auto C = m_dyn->Cdiag() * (1.0 / (theta * dt * dt));
    auto A = m_dyn->Adiag() * (1.0 / (theta * dt));

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    TripletMatrix P = TripletMatrix::fromDiag(diag);
    m_P_cache = P.asSparse();
    return m_P_cache;
}

VectorXr& ConicxxAssembler::q(real_t dt, real_t theta) {
    (void)theta;

    auto q_top = -(m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal()));
    auto q_bot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    m_q_cache.resize(q_top.size() + q_bot.size());
    m_q_cache << q_top, q_bot;

    return m_q_cache;
}

const SparseMatrix<Eigen::ColMajor>& ConicxxAssembler::A(real_t dt, real_t theta) {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * (1.0 / (theta * dt * dt)));
    TripletMatrix A_diag = TripletMatrix::fromDiag(m_dyn->Adiag() * (1.0 / (theta * dt)));

    const int n_contact = m_dyn->numContactRows();
    const int n_var = m_dyn->numV() + m_dyn->numSprings() + m_dyn->numDampers();
    VectorXr Smu = computeSmu();

    TripletMatrix A = (m_dyn->Wg() | (C * -1.0) | TripletMatrix::zero(C.nRows(), A_diag.nCols()))
                          .vConcat(m_dyn->Wgamma() | TripletMatrix::zero(A_diag.nRows(), C.nCols()) | (A_diag * -1.0))
                          .vConcat(m_dyn->W().scaleRows(Smu) * -1.0 | TripletMatrix::zero(n_contact, m_dyn->numSprings() + m_dyn->numDampers()));

    (void)n_var;
    m_A_cache = A.asSparse();
    return m_A_cache;
}

VectorXr& ConicxxAssembler::b(real_t dt, real_t theta) {
    auto lambda_g = m_dyn->Lambda_g();
    if (lambda_g.size() != m_dyn->Cdiag().size()) lambda_g = VectorXr::Zero(m_dyn->Cdiag().size());

    const auto& Wcontact = m_dyn->W().asSparse();
    VectorXr Smu = computeSmu();
    VectorXr restitution = m_dyn->restitutionVec();
    VectorXr biasImpulse = m_dyn->contactVVec();

    VectorXr gamma_old = Wcontact * m_dyn->vVec() + m_dyn->contactVVec();

    if ((restitution.array() != 0.0).any())
        biasImpulse += restitution.cwiseProduct(gamma_old);

    VectorXr b_top = -(1.0 / (theta * dt * dt)) * m_dyn->Cdiag().cwiseProduct(lambda_g) - ((1.0 - theta) / theta) * (m_dyn->Wg().asSparse() * m_dyn->vVec()) - (1.0 / theta) * m_dyn->C_v_vec();
    VectorXr b_mid = -((1.0 - theta) / theta) * (m_dyn->Wgamma().asSparse() * m_dyn->vVec()) - (1.0 / theta) * m_dyn->A_v_vec();
    VectorXr b_contact = Smu.cwiseProduct(biasImpulse);

    // De Saxce shift s_i = mu_i * ||u_T(v^n)|| on the (unscaled) normal row.
    // Under P_mu = diag(1, mu, mu) the normal row carries physical units, so the
    // shift keeps its physical factor mu.
    for (int i = m_dyn->numFrictionlessContacts(); i < m_dyn->numContactRows(); i += 3) {
        const real_t y_1 = gamma_old[i + 1];
        const real_t y_2 = gamma_old[i + 2];
        b_contact[i] += m_dyn->muVec()[i] * std::sqrt(y_1 * y_1 + y_2 * y_2);
    }

    auto beta = m_dyn->system().config().constraint_bias_factor;
    if (beta > 0) b_top.noalias() -= (-m_dyn->Cdiag().cwiseProduct(lambda_g) / dt + m_dyn->g_error_vec()) * (beta / (dt * theta));

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

VectorXr ConicxxAssembler::computeSmu() const {
    VectorXr Smu = VectorXr::Ones(m_dyn->numContactRows());
    for (const auto& c : m_dyn->contacts()) {
        // friction less rows
        if (c.friction_mu <= 0) continue;
        if (c.impulse_size != 3) continue;

        // // old convention with 1 / mu
        // const int idx = c.impulse_base_index;
        // if (idx < 0 || idx >= Smu.size()) continue;
        // Smu[idx] = 1.0 / c.friction_mu;

        // Scaling convention P_mu = diag(1, mu, mu) per frictional contact
        // (Acary et al. 2024, Sec. 3): tangential rows are scaled DOWN by mu,
        // the normal row stays at 1. Physical impulse recovery: lambda = Smu * z.
        const int idx = c.impulse_base_index;
        if (idx < 0 || idx + 2 >= Smu.size()) continue;
        Smu[idx + 1] = c.friction_mu;
        Smu[idx + 2] = c.friction_mu;
    }
    return Smu;
}

}  // namespace cardillo::physics::assembly
