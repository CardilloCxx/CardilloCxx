#include "pgs_assembler.hpp"

#include <Eigen/Cholesky>
#include <cmath>
#include <optional>

#include "../../rigid_body/rigid_body.hpp"
#include "../constraints/constraints.hpp"

namespace cardillo::physics::assembly {

VectorXr PgsAssembler::rhs(real_t dt, real_t theta, const VectorXr& u_free) const {
    const auto& vn = m_dyn->vVec();
    const auto& Wg = m_dyn->Wg().asSparse();
    const auto& Wgamma = m_dyn->Wgamma().asSparse();
    const auto& Wcontact = m_dyn->W().asSparse();
    const auto& M_diag = m_dyn->MDiag();
    const auto& M_inv = m_dyn->MinvDiag();
    const auto& restitution = m_dyn->restitutionVec();

    const int nSprings = m_dyn->numSprings();
    const int nDampers = m_dyn->numDampers();
    const int nContacts = m_dyn->numContactRows();
    const int nLambda = nSprings + nDampers + nContacts;

    VectorXr rhs_vel = M_diag.cwiseProduct(vn) + dt * m_dyn->fVecExternal();
    if (!m_cfg.moreau_implicit_gyroscopy)
        rhs_vel += dt * m_dyn->fVecGyroscopic();
    else
        std::cerr << "PgsAssembler::rhs: Warning: Gyroscopic forces are treated implicitly, not implemented in PGS assembler rhs\n";

    VectorXr rhs = VectorXr::Zero(nLambda);
    if (nSprings > 0) {
        rhs.head(nSprings) = m_dyn->springBias(dt, theta);
        rhs.head(nSprings).noalias() += Wg * M_inv.cwiseProduct(rhs_vel);
    }

    if (nDampers > 0) {
        rhs.segment(nSprings, nDampers) = m_dyn->damperBias(dt, theta);

        rhs.segment(nSprings, nDampers).noalias() += Wgamma * M_inv.cwiseProduct(rhs_vel);
    }

    if (nContacts > 0) {
        VectorXr biasImpulse = m_dyn->contactVVec();

        if ((restitution.array() != 0.0).any()) {
            VectorXr v_prev = Wcontact * m_dyn->vVec() + m_dyn->contactVVec();
            biasImpulse += restitution.cwiseProduct(v_prev);
        }
        
        rhs.segment(nSprings + nDampers, nContacts).noalias() = Wcontact * u_free + biasImpulse; 
    }

    return rhs;
}

VectorXr PgsAssembler::ufree(real_t dt, real_t /*theta*/) const {
    const auto& vn = m_dyn->vVec();
    const auto& M_inv = m_dyn->MinvDiag();

    VectorXr ufree = vn + M_inv.cwiseProduct(dt * m_dyn->fVecExternal());
    if (!m_cfg.moreau_implicit_gyroscopy) ufree += M_inv.cwiseProduct(dt * m_dyn->fVecGyroscopic());
    return ufree;
}

Eigen::MatrixXd selectCols(const Eigen::MatrixXd& A, const std::vector<bool>& mask) {
    int rows = A.rows();
    int cols = A.cols();

    assert((int)mask.size() == cols && "Mask size must match the number of columns in A");

    // count selected cols
    int newCols = 0;
    for (bool b : mask)
        if (b) ++newCols;

    Eigen::MatrixXd B(rows, newCols);

    int j = 0;
    for (int i = 0; i < cols; ++i) {
        if (mask[i]) {
            B.col(j++) = A.col(i);
        }
    }

    return B;
}

BlockDiagonal PgsAssembler::DinvDiag(real_t dt, real_t theta) const {
    const int nSprings = m_dyn->numSprings();
    const int nDampers = m_dyn->numDampers();
    const int nContacts = m_dyn->numContactRows();
    const int nFrictionless = m_dyn->numFrictionlessContacts();

    auto diag_sparse = [](const Eigen::SparseMatrix<real_t, Eigen::RowMajor>& W, const VectorXr& MinvDiag) {
        const int C = (int)W.rows();
        VectorXr D = VectorXr::Zero(C);
        for (int cid = 0; cid < C; ++cid) {
            real_t Dii = 0;
            for (Eigen::SparseMatrix<real_t, Eigen::RowMajor>::InnerIterator it(W, cid); it; ++it) {
                const real_t w = it.value();
                Dii += w * w * MinvDiag[it.col()];
            }
            D[cid] = Dii;
        }
        return D;
    };

    VectorXr D_diag = VectorXr::Zero(nSprings + nDampers + nContacts);
    D_diag.head(nSprings) = diag_sparse(m_dyn->Wg().asSparse(), m_dyn->MinvDiag());
    D_diag.head(nSprings) += m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta);

    D_diag.segment(nSprings, nDampers) = diag_sparse(m_dyn->Wgamma().asSparse(), m_dyn->MinvDiag());
    D_diag.segment(nSprings, nDampers) += m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta);

    D_diag.segment(nSprings + nDampers, nContacts) = diag_sparse(m_dyn->W().asSparse(), m_dyn->MinvDiag());

    BlockDiagonal Dinv;
    int offset = 0;
    const auto& reg = m_dyn->system().ecs();

    for (const auto& constraint : m_dyn->constraintResults()) {
        const bool addA = !RigidBody::isStatic(reg, constraint.a);
        const bool addB = !RigidBody::isStatic(reg, constraint.b);
        if (!addA && !addB) continue;

        const auto& c_used = constraint.c_used;
        const int nSp = std::count(c_used.begin(), c_used.end(), true);

        if (nSp > 0) Dinv.addBlockDiag(D_diag.segment(offset, nSp).cwiseInverse());
        offset += nSp;
    }

    for (const auto& constraint : m_dyn->constraintResults()) {
        const bool addA = !RigidBody::isStatic(reg, constraint.a);
        const bool addB = !RigidBody::isStatic(reg, constraint.b);
        if (!addA && !addB) continue;

        const auto& a_used = constraint.a_used;
        const int nDa = std::count(a_used.begin(), a_used.end(), true);

        if (nDa > 0) Dinv.addBlockDiag(D_diag.segment(offset, nDa).cwiseInverse());
        offset += nDa;
    }

    for (int cid = 0; cid < nFrictionless; ++cid) {
        real_t Dii = D_diag[offset + cid];
        Dinv.addBlockDiag(VectorXr::Constant(1, Dii > 0 ? 1.0 / Dii : 0));
    }

    for (int cid = nFrictionless; cid < nFrictionless + 3 * m_dyn->numFrictionalContacts(); cid += 3) {
        VectorXr Dii = D_diag.segment(offset + cid, 3);
        Dinv.addBlockDiag(Dii.cwiseInverse());
    }
    Dinv.calcDiag();
    return Dinv;
}

Eigen::SparseMatrix<real_t> PgsAssembler::W() const {
    return m_dyn->Wg().vConcat(m_dyn->Wgamma()).vConcat(m_dyn->W()).asSparse();
}

VectorXr PgsAssembler::C(real_t dt, real_t theta) const {
    const int nSprings = m_dyn->numSprings();
    const int nDampers = m_dyn->numDampers();
    const int nContacts = m_dyn->numContactRows();
    const int nLambda = nSprings + nDampers + nContacts;

    VectorXr C_vec = VectorXr::Zero(nLambda);
    if (nSprings > 0) C_vec.head(nSprings) = m_dyn->Cdiag() * DynamicsAssembler::springComplianceScale(dt, theta);
    if (nDampers > 0) C_vec.segment(nSprings, nDampers) = m_dyn->Adiag() * DynamicsAssembler::damperComplianceScale(dt, theta);
    if (nContacts > 0) C_vec.tail(nContacts) = VectorXr::Zero(nContacts);

    return C_vec;
}

}  // namespace cardillo::physics::assembly