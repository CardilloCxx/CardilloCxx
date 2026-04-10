#include "qoco_assembler.hpp"
#include <stdexcept>
#include "dynamics_assembler.hpp"

namespace cardillo::physics::assembly {

QOCOCscMatrix* QocoAssembler::toQocoCSC(SparseMatrix<Eigen::ColMajor>& A, QOCOCscMatrix& out) {
    using StorageIndex = typename SparseMatrix<Eigen::ColMajor>::StorageIndex;
    static_assert(sizeof(StorageIndex) == sizeof(QOCOInt),
                  "SparseMatrix<>::StorageIndex and QOCOInt must have the same size for zero-copy "
                  "conversion");
    static_assert(sizeof(typename SparseMatrix<Eigen::ColMajor>::Scalar) == sizeof(QOCOFloat), "SparseMatrix<>::Scalar and QOCOFloat must have the same size for zero-copy conversion");
    out.m = static_cast<QOCOInt>(A.rows());
    out.n = static_cast<QOCOInt>(A.cols());
    out.nnz = static_cast<QOCOInt>(A.nonZeros());
    out.p = const_cast<QOCOInt*>(A.outerIndexPtr());
    out.i = const_cast<QOCOInt*>(A.innerIndexPtr());
    out.x = const_cast<QOCOFloat*>(A.valuePtr());

    return &out;
}

QOCOFloat* QocoAssembler::toQocoVector(Eigen::VectorX<real_t>& v) {
    static_assert(sizeof(real_t) == sizeof(QOCOFloat), "real_t and QOCOFloat must have the same size for zero-copy conversion");
    return reinterpret_cast<QOCOFloat*>(v.data());
}

QOCOCscMatrix* QocoAssembler::P(real_t dt, real_t theta) {
    auto M = m_dyn->MDiag();
    auto C = m_dyn->Cdiag() * (1.0 / (theta * dt * dt));
    auto A = m_dyn->Adiag() * (1.0 / (theta * dt));

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    TripletMatrix P = TripletMatrix::fromDiag(diag);
    m_P_cache = P.asSparse();

    return toQocoCSC(m_P_cache, m_P_view);
}

QOCOCscMatrix* QocoAssembler::A(real_t dt, real_t theta) {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * (1.0 / (theta * dt * dt)));
    TripletMatrix A = TripletMatrix::fromDiag(m_dyn->Adiag() * (1.0 / (theta * dt)));

    A = (m_dyn->Wg() | (C * -1.0) | TripletMatrix::zero(C.nRows(), A.nCols())).vConcat(m_dyn->Wgamma() | TripletMatrix::zero(A.nRows(), C.nCols()) | (A * -1.0));

    m_A_cache = A.asSparse();
    return toQocoCSC(m_A_cache, m_A_view);
}

VectorXr QocoAssembler::computeSmu() {
    VectorXr Smu = VectorXr::Ones(m_dyn->numContactRows());
    for (const auto& c : m_dyn->contacts()) {
        if (c.friction_mu <= 0) continue;
        const int idx = c.impulse_base_index;
        if (idx < 0 || idx >= Smu.size()) continue;
        Smu[idx] = 1.0 / c.friction_mu;
    }
    return Smu;
}

QOCOCscMatrix* QocoAssembler::G(real_t dt, real_t theta) {
    VectorXr Smu = computeSmu();
    TripletMatrix G = ((m_dyn->W().scaleRows(Smu) * -1.0) | TripletMatrix::zero(m_dyn->numContactRows(), m_dyn->numSprings() + m_dyn->numDampers()));

    m_G_cache = G.asSparse();
    return toQocoCSC(m_G_cache, m_G_view);
}

QOCOFloat* QocoAssembler::c(real_t dt, real_t theta) {
    auto cTop = -(m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal()));
    auto cBot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    m_c_cache.resize(cTop.size() + cBot.size());
    m_c_cache << cTop, cBot;
    return toQocoVector(m_c_cache);
}

QOCOFloat* QocoAssembler::b(real_t dt, real_t theta) {
    auto Lambda_g = m_dyn->Lambda_g();
    if (Lambda_g.size() != m_dyn->Cdiag().size()) {
        Lambda_g = VectorXr::Zero(m_dyn->Cdiag().size());
    }

    auto bTop = -(1.0 / (theta * dt * dt)) * m_dyn->Cdiag().cwiseProduct(Lambda_g) - ((1.0 - theta) / theta) * m_dyn->Wg().asSparse() * m_dyn->vVec();
    auto bBot = -((1.0 - theta) / theta) * m_dyn->Wgamma().asSparse() * m_dyn->vVec();

    m_b_cache.resize(bTop.size() + bBot.size());
    m_b_cache << bTop, bBot;
    return toQocoVector(m_b_cache);
}

QOCOFloat* QocoAssembler::h(real_t dt, real_t theta) {
    // Smu missing, for w = 0
    const int n = m_dyn->numContactRows();
    m_h_cache = VectorXr::Zero(n);
    return toQocoVector(m_h_cache);
}

}  // namespace cardillo::physics::assembly
