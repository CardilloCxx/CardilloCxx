#include "qoco_assembler.hpp"
#include "dynamics_assembler.hpp"
#include <cstdlib>
#include <stdexcept>

namespace cardillo::physics::assembly {

QOCOCscMatrix QocoAssembler::toQocoCSC(const SparseMatrix<Eigen::ColMajor>& A) {
    QOCOCscMatrix out;

    out.m = A.rows();
    out.n = A.cols();
    out.nnz = A.nonZeros();

    out.p = (QOCOInt*)std::malloc(sizeof(QOCOInt) * (A.cols() + 1));
    out.i = (QOCOInt*)std::malloc(sizeof(QOCOInt) * A.nonZeros());
    out.x = (QOCOFloat*)std::malloc(sizeof(QOCOFloat) * A.nonZeros());

    std::memcpy(out.p, A.outerIndexPtr(), sizeof(QOCOInt) * (A.cols() + 1));
    std::memcpy(out.i, A.innerIndexPtr(), sizeof(QOCOInt) * A.nonZeros());
    std::memcpy(out.x, A.valuePtr(),      sizeof(QOCOFloat) * A.nonZeros());
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

QOCOVectorf* QocoAssembler::toQocoVectorf(const Eigen::VectorX<real_t>& v) {
    const int n = static_cast<int>(v.size());
    return new_qoco_vectorf(reinterpret_cast<const QOCOFloat*>(v.data()), n);
}

QOCOCscMatrix  QocoAssembler::P() {
    auto M = m_dyn->MDiag(); 
    auto C = m_dyn->Cdiag();
    auto A = m_dyn->Adiag();

    Eigen::VectorX<real_t> diag(M.size() + C.size() + A.size());
    diag << M, C, A;

    return toQocoCSC(diag);
}

QOCOCscMatrix QocoAssembler::A() {
    TripletMatrix C = TripletMatrix::fromDiag(m_dyn->Cdiag());
    TripletMatrix A = TripletMatrix::fromDiag(m_dyn->Adiag());

    A = (m_dyn->Wg().T()    | (C * -1.0)   | A.zeroLike()).vConcat(
        m_dyn->Wgamma().T() | C.zeroLike() | (A * -1.0));

    return toQocoCSC( A.asSparse() );
}

QOCOCscMatrix QocoAssembler::G() {
    // Smu missing
    TripletMatrix G = (m_dyn->W().T() | m_dyn->Wg().T().zeroLike() | m_dyn->Wgamma().T().zeroLike());

    return toQocoCSC( G.asSparse() );
}

QOCOVectorf* QocoAssembler::c() {
    auto cTop = - (m_dyn->MDiag().cwiseProduct(m_dyn->vVec()) + m_dyn->fVecGyroscopic());
    auto cBot = VectorXr::Zero(m_dyn->numSprings() + m_dyn->numDampers());

    Eigen::VectorX<real_t> c(cTop.size() + cBot.size());
    c << cTop, cBot;
    return toQocoVectorf(c);
}

QOCOVectorf* QocoAssembler::b() {

    auto bTop = -2.0 * m_dyn->Cdiag().cwiseProduct(m_dyn->Lambda_g()) - m_dyn->Wg().asSparse().transpose() * m_dyn->vVec();
    auto bBot = - m_dyn->Wgamma().asSparse().transpose() * m_dyn->vVec();

    Eigen::VectorX<real_t> b(bTop.size() + bBot.size());
    b << bTop, bBot;
    return toQocoVectorf(b);
}

QOCOVectorf* QocoAssembler::h() {
    //Smu missing, w always = 0
    const int n = m_dyn->W().nCols();
    auto h = VectorXr::Zero(n);
    return toQocoVectorf(h);
}

} // namespace cardillo::physics::assembly

