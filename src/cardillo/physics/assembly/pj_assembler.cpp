#include "pj_assembler.hpp"

#include <iostream>

namespace cardillo::physics::assembly {
PjAssembler::~PjAssembler() = default;
bool PjAssembler::buildAndFactorS(real_t dt, real_t theta, bool implicitGyro) {
    auto sc = m_dyn->timings()->scope(misc::TimingManager::TimerId::BuildAndFactorS);

    const int totalV = (m_dyn->bodyVelOffsets().empty() ? 0 : m_dyn->bodyVelOffsets().back());
    const CscMatrix& Wg = m_dyn->Wg();
    const CscMatrix& Wgamma = m_dyn->Wgamma();
    const int nSprings = (int)Wg.rows();
    const int nDampers = (int)Wgamma.rows();
    const int extV = totalV + nSprings + nDampers;

    if (extV == 0) {
        m_S_sparse_lu.reset();
        return true;
    }

    // Saddle-point matrix of the theory manual (eq:saddle_point) in (u, pi_g, pi_gamma), pi = h*lambda:
    //   [ M_eff  -W_g          -W_gamma        ]
    //   [ -W_g^T -C/(th^2 h^2)  0              ]
    //   [ -W_gamma^T  0        -A/(th h)       ]
    // The entries are added in the same order every step, zeros included; the pattern of S, and with
    // it the analysis of its factorization, changes only when an entry becomes nonzero for the first
    // time or the structure changes (see misc::FixedPatternMatrix).
    m_S.begin(extV);
    for (int i = 0; i < totalV; ++i) m_S.add(i, i, m_dyn->MDiag()[i]);

    // Optionally include linearized gyroscopic term in the system matrix.
    // This adds -dt * G(omega_n) to each rigid body's rotational 3x3 block and breaks symmetry.
    if (implicitGyro) {
        const auto& reg = m_dyn->system().ecs();
        auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b < 0 || b + 1 >= (int)m_dyn->bodyVelOffsets().size()) continue;
            const int off = m_dyn->bodyVelOffsets()[(size_t)b];
            const int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - off;
            if (nV < 6) continue;
            if (!reg.all_of<C_RigidBodyTag, C_AngularVelocity3>(e)) continue;
            const Vector3r& omega = reg.get<C_AngularVelocity3>(e).value;
            const Vector3r& I = m_dyn->system().getInertiaDiag(e);
            const Vector3r& Iomega = I.cwiseProduct(omega);
            const Matrix33r omegaSkew = SkewSymmetricMatrix3r(omega);
            const Matrix33r IomegaSkew = SkewSymmetricMatrix3r(Iomega);

            // G_rot = 0.5 * ( [I*omega]_x - [omega]_x * I )
            const Matrix33r Grot = (real_t)0.5 * (IomegaSkew - omegaSkew * I.asDiagonal());
            const Matrix33r corr = -dt * Grot;  // M_eff = M - dt*G
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) m_S.add(off + 3 + r, off + 3 + c, corr(r, c));
        }
    }

    // -W_g / -W_gamma below and their transposes beside the mass block
    auto addRows = [&](const CscMatrix& rows, int rowOffset) {
        for (int col = 0; col < rows.outerSize(); ++col) {
            for (CscMatrix::InnerIterator it(rows, col); it; ++it) {
                m_S.add(rowOffset + it.row(), col, -it.value());
                m_S.add(col, rowOffset + it.row(), -it.value());
            }
        }
    };
    addRows(Wg, totalV);
    addRows(Wgamma, totalV + nSprings);
    const real_t cScale = -DynamicsAssembler::springComplianceScale(dt, theta);
    const real_t aScale = -DynamicsAssembler::damperComplianceScale(dt, theta);
    for (int i = 0; i < nSprings; ++i) m_S.add(totalV + i, totalV + i, cScale * m_dyn->Cdiag()[i]);
    for (int i = 0; i < nDampers; ++i) m_S.add(totalV + nSprings + i, totalV + nSprings + i, aScale * m_dyn->Adiag()[i]);
    const bool patternChanged = m_S.finish(extV);
    const CscMatrix& S = m_S.matrix();

    // Factorize using SparseLU; analyzePattern() (fill-reducing ordering, symbolic elimination) only
    // after a change of the pattern.
    const bool wantSymmetric = !implicitGyro;
    const bool analyze = patternChanged || !m_S_sparse_lu.has_value() || wantSymmetric != m_analyzedSymmetric;
    if (analyze && m_cfg.solver_verbosity >= 2) std::cout << "[PJ] analyzing the pattern of S: n = " << S.rows() << ", nnz = " << S.nonZeros() << std::endl;
    try {
        if (analyze) {
            m_S_sparse_lu.emplace();
            m_S_sparse_lu->isSymmetric(wantSymmetric);
            m_S_sparse_lu->analyzePattern(S);
            m_analyzedSymmetric = wantSymmetric;
        }
        m_S_sparse_lu->factorize(S);
        if (m_S_sparse_lu->info() != Eigen::Success) {
            m_S_sparse_lu.reset();
            std::cout << "PjAssembler::buildAndFactorS: SparseLU factorization failed\n";
            return false;
        } else if (m_cfg.debug_rb) {
            std::cout << "[PjAssembler] SparseLU factorization success.\n";
        }
    } catch (const std::exception& ex) {
        if (m_cfg.debug_rb) {
            std::cout << "[PjAssembler] Exception during SparseLU: " << ex.what() << '\n';
        }
        m_S_sparse_lu.reset();
        return false;
    }
    return true;
}

// RHS
VectorXr PjAssembler::rhs(real_t dt, real_t theta) const {
    const auto& vn = m_dyn->vVec();
    const auto& fn_ext = m_dyn->fVecExternal();     // gravity + applied external forces
    const auto& fn_gyro = m_dyn->fVecGyroscopic();  // gyroscopic forces from current state
    const bool implicitGyro = m_cfg.moreau_implicit_gyroscopy;
    const auto& M_diag = m_dyn->MDiag();
    const int totalV = (m_dyn->bodyVelOffsets().empty() ? 0 : m_dyn->bodyVelOffsets().back());
    const int nSprings = (int)m_dyn->Cdiag().size();
    const int nDampers = (int)m_dyn->Adiag().size();
    const int extV = totalV + nSprings + nDampers;

    // RHS: M*vn + dt*f_ext (+ dt*f_gyro if treated explicitly)
    VectorXr rhs = VectorXr::Zero((index_t)extV);
    rhs.segment(0, totalV) = M_diag.cwiseProduct(vn) + dt * fn_ext;
    if (!implicitGyro) rhs.segment(0, totalV) += dt * fn_gyro;
    if (nSprings > 0) rhs.segment(totalV, nSprings) = m_dyn->springBias(dt, theta);
    if (nDampers > 0) rhs.segment(totalV + nSprings, nDampers) = m_dyn->damperBias(dt, theta);

    return rhs;
}

// Solve full extended system and return complete solution
VectorXr PjAssembler::solveS(const VectorXr& rhs_ext) const {
    if (rhs_ext.size() == 0) return rhs_ext;  // empty system: identity solve
    if (!m_S_sparse_lu.has_value()) throw std::runtime_error("PjAssembler::solveS called but S matrix is not factorized");

    Eigen::VectorXd sol = m_S_sparse_lu->solve(rhs_ext);
    if (m_S_sparse_lu->info() != Eigen::Success) {
        throw std::runtime_error("PjAssembler::solveS: SparseLU solve failed");
    }
    return sol;
}

MatrixXXr PjAssembler::solveSColumns(const MatrixXXr& B) const {
    if (B.size() == 0) return B;
    if (!m_S_sparse_lu.has_value()) throw std::runtime_error("PjAssembler::solveSColumns called but S matrix is not factorized");

    MatrixXXr X = m_S_sparse_lu->solve(B);
    if (m_S_sparse_lu->info() != Eigen::Success) {
        throw std::runtime_error("PjAssembler::solveSColumns: SparseLU solve failed");
    }
    return X;
}
}  // namespace cardillo::physics::assembly