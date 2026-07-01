#include "pj_assembler.hpp"

namespace cardillo::physics::assembly {
PjAssembler::~PjAssembler() = default;
bool PjAssembler::buildAndFactorS(real_t dt, real_t theta, bool implicitGyro, bool lambdaTheta) {
    auto sc = m_dyn->timings()->scope(cardillo::misc::TimingManager::TimerId::BuildAndFactorS);

    const int totalV = (m_dyn->bodyVelOffsets().empty() ? 0 : m_dyn->bodyVelOffsets().back());
    const int nSprings = m_dyn->Wg().nRows();
    const int nDampers = m_dyn->Wgamma().nRows();
    const int extV = totalV + nSprings + nDampers;

    if (extV == 0) {
        m_S_sparse_lu.reset();
        m_S = TripletMatrix::zero(0, 0);
        return true;
    }

    std::vector<Eigen::Triplet<real_t>> mTrips;
    mTrips.reserve(static_cast<std::size_t>(totalV) + 64);

    // Top-left: M diagonal
    for (int i = 0; i < totalV; ++i) {
        real_t mval = m_dyn->MDiag()[i];
        if (mval != (real_t)0) mTrips.emplace_back(i, i, mval);
    }

    // Optionally include linearized gyroscopic term in the system matrix.
    // This adds -dt * G(omega_n) to each rigid body's rotational 3x3 block and breaks symmetry.
    if (implicitGyro) {
        const auto& reg = m_dyn->system().ecs();
        auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b < 0 || b + 1 >= (int)m_dyn->bodyVelOffsets().size()) continue;
            const int off = m_dyn->bodyVelOffsets()[(size_t)b];
            const int nV = m_dyn->bodyVelOffsets()[(size_t)b + 1] - off;
            if (nV < 6) continue;
            if (!reg.all_of<cardillo::C_RigidBodyTag, cardillo::C_AngularVelocity3>(e)) continue;
            const Vector3r omega = reg.get<cardillo::C_AngularVelocity3>(e).value;  // body-frame
            const Vector3r I = m_dyn->system().getInertiaDiag(e);
            const Vector3r Iomega = I.cwiseProduct(omega);
            const Matrix33r Idiag = I.asDiagonal().toDenseMatrix();
            const Matrix33r omegaSkew = skew_from_vector(omega);
            const Matrix33r IomegaSkew = skew_from_vector(Iomega);

            // G_rot = 0.5 * ( [I*omega]_x - [omega]_x * I )
            const Matrix33r Grot = (real_t)0.5 * (IomegaSkew - omegaSkew * Idiag);
            const Matrix33r corr = -dt * Grot;  // M_eff = M - dt*G

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    const real_t val = corr(r, c);
                    if (val != (real_t)0) mTrips.emplace_back(off + 3 + r, off + 3 + c, val);
                }
            }
        }
    }

    const real_t topScale = lambdaTheta ? theta : (real_t)1.0;

    const TripletMatrix Mblk(totalV, totalV, std::make_shared<std::vector<Eigen::Triplet<real_t>>>(std::move(mTrips)));
    const real_t cScale = -(real_t)1.0 / (theta * dt * dt);
    const real_t aScale = -(real_t)1.0 / (theta * dt);

    const TripletMatrix top = Mblk | (m_dyn->Wg() * topScale).T() | (m_dyn->Wgamma() * topScale).T();
    const TripletMatrix mid = m_dyn->Wg() | (TripletMatrix::fromDiag(m_dyn->Cdiag()) * cScale) | TripletMatrix::zero(nSprings, nDampers);
    const TripletMatrix bot = m_dyn->Wgamma() | TripletMatrix::zero(nDampers, nSprings) | (TripletMatrix::fromDiag(m_dyn->Adiag()) * aScale);

    m_S = top.vConcat(mid).vConcat(bot);
    const auto& S_sparse = m_S.asSparse();

    // Factorize using SparseLU
    try {
        m_S_sparse_lu.emplace();
        m_S_sparse_lu->isSymmetric(!implicitGyro && !lambdaTheta);
        m_S_sparse_lu->analyzePattern(S_sparse);
        m_S_sparse_lu->factorize(S_sparse);
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
    const bool lambdaTheta = m_cfg.moreau_lambda_theta;
    const auto& Wg = m_dyn->Wg().asSparse();
    const auto& Wgamma = m_dyn->Wgamma().asSparse();
    const auto& M_diag = m_dyn->MDiag();
    const int totalV = (m_dyn->bodyVelOffsets().empty() ? 0 : m_dyn->bodyVelOffsets().back());
    const int nSprings = (int)m_dyn->Cdiag().size();
    const int nDampers = (int)m_dyn->Adiag().size();
    const int extV = totalV + nSprings + nDampers;
    const auto& Cdiag = m_dyn->Cdiag();
    const auto& C_v_vec = m_dyn->C_v_vec();
    const auto& A_v_vec = m_dyn->A_v_vec();

    // Lambda vectors may be uninitialized on first step; copy locally and ensure correct sizes
    VectorXr Lambda_g = m_dyn->Lambda_g();
    if ((int)Lambda_g.size() != nSprings) Lambda_g = VectorXr::Zero(nSprings);
    VectorXr Lambda_gamma = m_dyn->Lambda_gamma();
    if ((int)Lambda_gamma.size() != nDampers) Lambda_gamma = VectorXr::Zero(nDampers);

    // RHS: M*vn + dt*f_ext (+ dt*f_gyro if treated explicitly)
    VectorXr rhs = VectorXr::Zero((index_t)extV);
    rhs.segment(0, totalV) = M_diag.cwiseProduct(vn) + dt * fn_ext;
    if (!implicitGyro) rhs.segment(0, totalV) += dt * fn_gyro;
    if (lambdaTheta && (nSprings > 0 || nDampers > 0)) {
        VectorXr corr = VectorXr::Zero(totalV);
        if (nSprings > 0) corr.noalias() += Wg.transpose() * Lambda_g;
        if (nDampers > 0) corr.noalias() += Wgamma.transpose() * Lambda_gamma;
        rhs.segment(0, totalV).noalias() -= (1.0 - theta) * corr;
    }
    if (nSprings > 0) rhs.segment(totalV, nSprings) = -(1.0 / (theta * dt * dt)) * Cdiag.cwiseProduct(Lambda_g) - ((1.0 - theta) / theta) * (Wg * vn) - (1.0 / theta) * C_v_vec;
    if (nDampers > 0) rhs.segment(totalV + nSprings, nDampers) = -((1.0 - theta) / theta) * (Wgamma * vn) - (1.0 / theta) * A_v_vec;

    auto beta = m_dyn->system().config().constraint_bias_factor;
    if (beta > 0) rhs.segment(totalV, nSprings).noalias() -= (-m_dyn->Cdiag().cwiseProduct(Lambda_g) / dt + m_dyn->g_error_vec()) * (beta / (dt * theta));

    // // TODO: Constraint stabilization proposed by Marco; this is not working, but I'm not sure why.
    // if (nSprings > 0) rhs.segment(totalV, nSprings) = (1.0 / (theta * dt * dt)) * m_dyn->g_error_vec() - ((1.0 - theta) / theta) * (Wg * vn) - (1.0 / theta) * C_v_vec;

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
}  // namespace cardillo::physics::assembly