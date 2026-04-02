#include "qoco_assembler.hpp"
#include "dynamics_assembler.hpp"
#include <cstdlib>
#include <stdexcept>

namespace cardillo::physics::assembly {

// QOCOCscMatrix QocoAssembler::toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A) {
//     QOCOCscMatrix out;
// 
//     out.m = A.rows();
//     out.n = A.cols();
//     out.nnz = A.nonZeros();
// 
//     out.p = (QOCOInt*)std::malloc(sizeof(QOCOInt) * (A.cols() + 1));
//     out.i = (QOCOInt*)std::malloc(sizeof(QOCOInt) * A.nonZeros());
//     out.x = (QOCOFloat*)std::malloc(sizeof(QOCOFloat) * A.nonZeros());
// 
//     std::memcpy(out.p, A.outerIndexPtr(), sizeof(QOCOInt) * (A.cols() + 1));
//     std::memcpy(out.i, A.innerIndexPtr(), sizeof(QOCOInt) * A.nonZeros());
//     std::memcpy(out.x, A.valuePtr(),      sizeof(QOCOFloat) * A.nonZeros());
//     return out;
// }

QOCOCscMatrix QocoAssembler::toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A) {
    QOCOCscMatrix out;

    out.m = static_cast<QOCOInt>(A.rows());
    out.n = static_cast<QOCOInt>(A.cols());
    out.nnz = static_cast<QOCOInt>(A.nonZeros());

    out.p = static_cast<QOCOInt*>(std::malloc(sizeof(QOCOInt) * (out.n + 1)));
    out.i = static_cast<QOCOInt*>(std::malloc(sizeof(QOCOInt) * out.nnz));
    out.x = static_cast<QOCOFloat*>(std::malloc(sizeof(QOCOFloat) * out.nnz));

    if (!out.p || !out.i || !out.x) {
        throw std::bad_alloc();
    }

    for (int k = 0; k < out.n + 1; ++k) {
        out.p[k] = static_cast<QOCOInt>(A.outerIndexPtr()[k]);
    }

    for (int k = 0; k < out.nnz; ++k) {
        out.i[k] = static_cast<QOCOInt>(A.innerIndexPtr()[k]);
    }

    for (int k = 0; k < out.nnz; ++k) {
        out.x[k] = static_cast<QOCOFloat>(A.valuePtr()[k]);
    }

    return out;
}

QOCOCscMatrix QocoAssembler::toQocoCSC(const Eigen::VectorX<real_t>& v) {
    const int n = static_cast<int>(v.size());

    QOCOCscMatrix out;
    out.m = out.n = n;
    out.nnz = n;

    out.p = (QOCOInt*)std::malloc(sizeof(QOCOInt) * (n + 1));
    out.i = (QOCOInt*)std::malloc(sizeof(QOCOInt) * n);
    out.x = (QOCOFloat*)std::malloc(sizeof(QOCOFloat) * n);

    if (!out.p || !out.i || !out.x) {
        throw std::bad_alloc();
    }

    for (int j = 0; j < n; ++j) {
        out.p[j] = j;       // one entry per column
        out.i[j] = j;       // diagonal
        out.x[j] = static_cast<QOCOFloat>(v[j]);
    }
    out.p[n] = n;

    return out;
}

QOCOFloat* QocoAssembler::toQocoVector(const Eigen::VectorX<real_t>& v) {
    const int n = static_cast<int>(v.size());
    QOCOFloat* out = static_cast<QOCOFloat*>(std::malloc(sizeof(QOCOFloat) * n));
    if (!out) {
        throw std::bad_alloc();
    }

    for (int i = 0; i < n; ++i) {
        out[i] = static_cast<QOCOFloat>(v[i]);
    }

    return out;
}

// C -> 4/dt^2 * C, 
// A -> 2/dt * A

QOCOCscMatrix  QocoAssembler::P(real_t dt, real_t theta) {
    auto M = m_dyn->MDiag(); 
    auto C = m_dyn->Cdiag() * (1.0 / (theta * dt * dt));
    auto A = m_dyn->Adiag() * (1.0 / (theta * dt));

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    return toQocoCSC(diag);
}

QOCOCscMatrix QocoAssembler::A(real_t dt, real_t theta) {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag() * (1.0 / (theta * dt * dt)));
    TripletMatrix A = TripletMatrix::fromDiag(m_dyn->Adiag() * (1.0 / (theta * dt)));

    A = (m_dyn->Wg()    | (C * -1.0)   | TripletMatrix::zero(C.nRows(), A.nCols())).vConcat(
        m_dyn->Wgamma() | TripletMatrix::zero(A.nRows(), C.nCols()) | (A * -1.0));

    return toQocoCSC( A.asSparse() );
}

QOCOCscMatrix QocoAssembler::G(real_t dt, real_t theta) {
    //TODO: add friction scaling, currently mu = 1.0
    TripletMatrix G = (m_dyn->W() | TripletMatrix::zero(m_dyn->W().nRows(), m_dyn->numSprings() + m_dyn->numDampers()));

    return toQocoCSC( G.asSparse() );
}

QOCOFloat* QocoAssembler::c(real_t dt, real_t theta) {
    auto cTop = - (m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + dt * (m_dyn->fVecGyroscopic() + m_dyn->fVecExternal()));
    auto cBot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    Eigen::VectorX<real_t> c(cTop.size() + cBot.size());
    c << cTop, cBot;
    return toQocoVector(c);
}

QOCOFloat* QocoAssembler::b(real_t dt, real_t theta) {

    auto Lambda_g = m_dyn->Lambda_g();
    if (Lambda_g.size() != m_dyn->Cdiag().size()) {
        Lambda_g = VectorXr::Zero(m_dyn->Cdiag().size());
    }
    
    auto bTop = - (1.0 / (theta * dt * dt)) * m_dyn->Cdiag().cwiseProduct(Lambda_g) - ((1.0 - theta) / theta) *  m_dyn->Wg().asSparse() * m_dyn->vVec();
    auto bBot = - ((1.0 - theta) / theta) * m_dyn->Wgamma().asSparse() * m_dyn->vVec();

    Eigen::VectorX<real_t> b(bTop.size() + bBot.size());
    b << bTop, bBot;
    return toQocoVector(b);
}

QOCOFloat* QocoAssembler::h(real_t dt, real_t theta) {
    //Smu missing, for w = 0
    const int n = m_dyn->W().nRows();
    auto h = VectorXr::Zero(n);
    return toQocoVector(h);
}

} // namespace cardillo::physics::assembly

