#include "pgs_assembler.hpp"

#include <Eigen/Cholesky>
#include <cmath>
#include <optional>

#include "../../rigid_body/rigid_body.hpp"
#include "../constraints/constraints.hpp"

namespace cardillo::physics::assembly {

VectorXr PgsAssembler::rhs(real_t dt, real_t theta) const {
    const auto& vn = m_dyn->vVec();
    const auto& Wg = m_dyn->Wg().asSparse();
    const auto& Wgamma = m_dyn->Wgamma().asSparse();
    const auto& M_diag = m_dyn->MDiag();
    const auto& M_inv = m_dyn->MinvDiag();
    const int totalV = m_dyn->numV();
    const int nSprings = m_dyn->numSprings();
    const int nDampers = m_dyn->numDampers();

    VectorXr Lambda_g = m_dyn->Lambda_g();
    if ((int)Lambda_g.size() != nSprings) Lambda_g = VectorXr::Zero(nSprings);
    VectorXr Lambda_gamma = m_dyn->Lambda_gamma();
    if ((int)Lambda_gamma.size() != nDampers) Lambda_gamma = VectorXr::Zero(nDampers);

    VectorXr rhs_vel = M_diag.cwiseProduct(vn) + dt * m_dyn->fVecExternal();
    if (!m_cfg.moreau_implicit_gyroscopy)
        rhs_vel += dt * m_dyn->fVecGyroscopic();
    else
        std::cerr << "PgsAssembler::rhs: Warning: Gyroscopic forces are treated implicitly, not implemented in PGS assembler rhs\n";

    if (m_cfg.moreau_lambda_theta && (nSprings > 0 || nDampers > 0)) {
        VectorXr corr = VectorXr::Zero(totalV);
        if (nSprings > 0) corr.noalias() += Wg.transpose() * Lambda_g;
        if (nDampers > 0) corr.noalias() += Wgamma.transpose() * Lambda_gamma;
        rhs_vel -= (1.0 - theta) * corr;
    }

    VectorXr rhs = VectorXr::Zero(nSprings + nDampers);
    if (nSprings > 0) {
        rhs.head(nSprings) = +(1.0 / (theta * dt * dt)) * m_dyn->Cdiag().cwiseProduct(Lambda_g) + ((1.0 - theta) / theta) * (Wg * vn) + (1.0 / theta) * m_dyn->C_v_vec();

        auto beta = m_dyn->system().config().constraint_bias_factor;
        if (beta > 0) rhs.head(nSprings).noalias() += (-m_dyn->Cdiag().cwiseProduct(Lambda_g) / dt + m_dyn->g_error_vec()) * (beta / (dt * theta));

        rhs.head(nSprings).noalias() += Wg * M_inv.cwiseProduct(rhs_vel);
    }

    if (nDampers > 0) {
        rhs.tail(nDampers) = +((1.0 - theta) / theta) * (Wgamma * vn) + (1.0 / theta) * m_dyn->A_v_vec();

        rhs.tail(nDampers).noalias() += Wgamma * M_inv.cwiseProduct(rhs_vel);
    }

    return rhs;
}

VectorXr PgsAssembler::ufree(real_t dt, real_t theta) const {
    const auto& vn = m_dyn->vVec();
    const auto& M_inv = m_dyn->MinvDiag();
    const int totalV = m_dyn->numV();

    VectorXr ufree = vn + M_inv.cwiseProduct(dt * m_dyn->fVecExternal());
    if (!m_cfg.moreau_implicit_gyroscopy) ufree += M_inv.cwiseProduct(dt * m_dyn->fVecGyroscopic());
    return ufree;
}

Eigen::MatrixXd selectRows(const Eigen::MatrixXd& A, const std::vector<bool>& mask) {
    int rows = A.rows();
    int cols = A.cols();

    assert((int)mask.size() == rows && "Mask size must match the number of rows in A");

    // count selected rows
    int newRows = 0;
    for (bool b : mask)
        if (b) ++newRows;

    Eigen::MatrixXd B(newRows, cols);

    int j = 0;
    for (int i = 0; i < rows; ++i) {
        if (mask[i]) {
            B.row(j++) = A.row(i);
        }
    }

    return B;
}

BlockDiagonal PgsAssembler::Dinv(real_t dt, real_t theta) const {
    std::vector<VectorXr> CrowsList;
    std::vector<VectorXr> ArowsList;

    std::vector<MatrixXXr> springBlocks;
    std::vector<MatrixXXr> damperBlocks;

    for (auto& constraint : m_dyn->constraintResults()) {
        const auto& Crows = constraint.Crows;
        const auto& Arows = constraint.Arows;
        const auto& c_used = constraint.c_used;
        const auto& a_used = constraint.a_used;
        const int nSprings = std::count(c_used.begin(), c_used.end(), true);
        const int nDampers = std::count(a_used.begin(), a_used.end(), true);
        const auto& reg = m_dyn->system().ecs();

        const bool addA = !RigidBody::isStatic(reg, constraint.a);
        const bool addB = !RigidBody::isStatic(reg, constraint.b);
        if (!addA && !addB) continue;

        if (nSprings > 0) {
            MatrixXXr blockSpring = MatrixXXr::Zero(nSprings, nSprings);

            if (addA) {
                int b = reg.get<cardillo::C_BodyIndex>(constraint.a).b;
                int row0 = m_dyn->bodyVelOffsets()[(size_t)b];
                int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - row0;
                const MatrixXXr WgA_sel = selectRows(constraint.WgA, c_used);
                blockSpring.noalias() += WgA_sel.transpose() * m_dyn->MinvDiag().segment(row0, nV).asDiagonal() * WgA_sel;
            }

            if (addB) {
                int b = reg.get<cardillo::C_BodyIndex>(constraint.b).b;
                int row0 = m_dyn->bodyVelOffsets()[(size_t)b];
                int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - row0;
                const MatrixXXr WgB_sel = selectRows(constraint.WgB, c_used);
                blockSpring.noalias() += WgB_sel.transpose() * m_dyn->MinvDiag().segment(row0, nV).asDiagonal() * WgB_sel;
            }

            for (int i = 0, j = 0; i < (int)c_used.size(); ++i) {
                if (!c_used[i]) continue;

                blockSpring(j, j) += m_dyn->Cdiag()[(size_t)i] / (theta * dt * dt);
                ++j;
            }

            springBlocks.push_back(blockSpring);
        }

        if (nDampers > 0) {
            MatrixXXr blockDamper = MatrixXXr::Zero(nDampers, nDampers);

            if (addA) {
                int b = reg.get<cardillo::C_BodyIndex>(constraint.a).b;
                int row0 = m_dyn->bodyVelOffsets()[(size_t)b];
                int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - row0;
                const MatrixXXr WgammaA_sel = selectRows(constraint.WgammaA, a_used);
                blockDamper.noalias() += WgammaA_sel.transpose() * m_dyn->MinvDiag().segment(row0, nV).asDiagonal() * WgammaA_sel;
            }

            if (addB) {
                int b = reg.get<cardillo::C_BodyIndex>(constraint.b).b;
                int row0 = m_dyn->bodyVelOffsets()[(size_t)b];
                int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - row0;
                const MatrixXXr WgammaB_sel = selectRows(constraint.WgammaB, a_used);
                blockDamper.noalias() += WgammaB_sel.transpose() * m_dyn->MinvDiag().segment(row0, nV).asDiagonal() * WgammaB_sel;
            }

            for (int i = 0, j = 0; i < (int)a_used.size(); ++i) {
                if (!a_used[i]) continue;

                blockDamper(j, j) += m_dyn->Adiag()[(size_t)i] / (theta * dt);
                ++j;
            }

            damperBlocks.push_back(blockDamper);
        }
    }

    BlockDiagonal D;
    for (MatrixXXr& block : springBlocks) D.addBlock(block);
    for (MatrixXXr& block : damperBlocks) D.addBlock(block);

    return D.calculateInverse();
}

}  // namespace cardillo::physics::assembly