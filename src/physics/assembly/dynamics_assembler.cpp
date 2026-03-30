
#include "dynamics_assembler.hpp"
#include "../constraints/constraints.hpp"
#include "../../collision/collision_coal.hpp"
#include <Eigen/Cholesky>
#include <cmath>
#include <iostream>


namespace cardillo::physics {

namespace {
// Build a 6-DoF W row for a rigid body side with sign s (+1 for B, -1 for A)
static inline MatrixXXr buildWRowRigid(const Vector3r& n_world,
                                       const Vector3r& r_body,
                                       const Vector3r& n_body,
                                       real_t s) {
    MatrixXXr w(1,6); w.setZero();
    // Translational contribution: s * n
    w(0,0) = s * n_world.x();
    w(0,1) = s * n_world.y();
    w(0,2) = s * n_world.z();
    // Rotational contribution in body frame: r_body x (s * n_body)
    const Vector3r t_body = r_body.cross(s * n_body);
    w(0,3) = t_body.x();
    w(0,4) = t_body.y();
    w(0,5) = t_body.z();
    return w;
}

// Build a 3-DoF W row for a point-mass side with sign s (+1 for B, -1 for A)
static inline MatrixXXr buildWRowPoint(const Vector3r& n_world, real_t s) {
    MatrixXXr w(1,3); w.setZero();
    w(0,0) = s * n_world.x();
    w(0,1) = s * n_world.y();
    w(0,2) = s * n_world.z();
    return w;
}

} // anonymous namespace

void DynamicsAssembler::updateContactsFromSystem() {
    if (!m_collision_mgr) throw std::runtime_error("DynamicsAssembler::updateContactsFromSystem: no CollisionCoal provided");
    if (m_world.consumeStructureDirty()) m_collision_mgr->rebuild();
    m_collision_mgr->applyTransforms();
        m_contacts_ptr = &m_collision_mgr->detectAll();
}

void DynamicsAssembler::setLambda_g(const VectorXr& lam) {
    m_Lambda_g = lam;
    static bool headerPrinted = false;
    static int step = 0;
}

// ---------- Cached API (block-based) ----------

const VectorXr& DynamicsAssembler::qVec() { return m_q_vec; }
const VectorXr& DynamicsAssembler::vVec() { return m_v_vec; }
const VectorXr& DynamicsAssembler::fVec() { return m_f_vec; }
const VectorXr& DynamicsAssembler::fVecExternal() { return m_f_vec_external; }
const VectorXr& DynamicsAssembler::fVecGyroscopic() { return m_f_vec_gyroscopic; }

// ---------- Rebuild helpers ----------

void DynamicsAssembler::rebuildMass_() {
    const int Nb = m_world.numBodies();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_Minv_diag = VectorXr::Zero(totalV);
    m_M_diag = VectorXr::Zero(totalV);

    const auto& reg = m_world.ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b; if (b < 0 || b >= Nb) continue;
        // Get inverse mass diagonal directly
        VectorXr MinvDiag = m_world.getMassInverseDiag(e);
        const int off = m_body_vel_offsets[(size_t)b];
        const int n = (int)MinvDiag.size();
        for (int i = 0; i < n; ++i) {
            const real_t inv = MinvDiag[i];
            m_Minv_diag[off + i] = inv;
            m_M_diag[off + i] = (inv > (real_t)0) ? (real_t)1 / inv : (real_t)0;
        }
    }
}

void DynamicsAssembler::rebuildForces_() {
    const int Nb = m_world.numBodies();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_f_vec = VectorXr::Zero(totalV);
    m_f_vec_external = VectorXr::Zero(totalV);
    m_f_vec_gyroscopic = VectorXr::Zero(totalV);
    m_v_compat.assign((size_t)Nb, VectorXr());
    const auto& reg = m_world.ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b; if (b >= 0 && b < Nb) {
            const VectorXr fb_ext = m_world.getForceExternal(e);
            const VectorXr fb_gyro = m_world.getForceGyroscopic(e);
            const VectorXr fb = fb_ext + fb_gyro;
            const int off = m_body_vel_offsets[(size_t)b];
            const int n = (int)fb.size();
            if (n > 0) {
                std::copy(fb.data(), fb.data() + n, m_f_vec.data() + off);
                std::copy(fb_ext.data(), fb_ext.data() + n, m_f_vec_external.data() + off);
                std::copy(fb_gyro.data(), fb_gyro.data() + n, m_f_vec_gyroscopic.data() + off);
            }
        }
    }
    // One-shot: clear external force/torque components after assembling forces
    auto &reg_mut = const_cast<entt::registry&>(m_world.ecs());
    auto viewF = reg_mut.view<cardillo::C_ExternalForce>();
    for (auto e : viewF) { reg_mut.get<cardillo::C_ExternalForce>(e).f.setZero(); }
    auto viewT = reg_mut.view<cardillo::C_ExternalTorque>();
    for (auto e : viewT) { reg_mut.get<cardillo::C_ExternalTorque>(e).tau.setZero(); }
}

void DynamicsAssembler::loadStateFromSystem() {
    const int Nb = m_world.numBodies();
    const int totalQ = (m_body_pos_offsets.empty() ? 0 : m_body_pos_offsets.back());
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_q_vec = VectorXr::Zero(totalQ);
    m_v_vec = VectorXr::Zero(totalV);
    m_v_compat.assign((size_t)Nb, VectorXr()); // deprecated compatibility only
    const auto& reg = m_world.ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < Nb) {
            const VectorXr qb = m_world.getPosition(e);
            const VectorXr vb = m_world.getVelocity(e);
            const int offQ = m_body_pos_offsets[(size_t)b];
            const int nQ = (int)qb.size();
            if (nQ > 0) std::copy(qb.data(), qb.data() + nQ, m_q_vec.data() + offQ);
            const int offV = m_body_vel_offsets[(size_t)b];
            const int nV = (int)vb.size();
            if (nV > 0) std::copy(vb.data(), vb.data() + nV, m_v_vec.data() + offV);
            m_v_compat[(size_t)b] = vb; // deprecated
        }
    }
}

void DynamicsAssembler::writePositionToSystem(const VectorXr& q) {
    const auto& reg = m_world.ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < (int)m_body_pos_offsets.size()-1) {
            const int offQ = m_body_pos_offsets[(size_t)b];
            const int nQ = m_body_pos_offsets[(size_t)b+1] - offQ;
            VectorXr qb = (nQ>0) ? q.segment(offQ, nQ) : VectorXr(0);
            if (qb.size() >= 3 && reg.any_of<cardillo::C_Position3>(e)) {
                const_cast<cardillo::C_Position3&>(reg.get<cardillo::C_Position3>(e)).value = qb.head<3>();
            }
            if (qb.size() >= 7 && reg.any_of<cardillo::C_Orientation>(e)) {
                Quaternion4r qn(qb.tail<4>());
                const Quaternion4r q_ref = reg.get<cardillo::C_Orientation>(e).value;
                const_cast<cardillo::C_Orientation&>(reg.get<cardillo::C_Orientation>(e)).value = MathHelper::alignQuaternionTo(qn, q_ref);
            }
        }
    }
    m_world.markStateDirty();
    m_world.markForcesDirty();
}

void DynamicsAssembler::writeVelocityToSystem(const VectorXr& v) {
    const auto& reg = m_world.ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < (int)m_body_vel_offsets.size()-1) {
            const int offV = m_body_vel_offsets[(size_t)b];
            const int nV = m_body_vel_offsets[(size_t)b+1] - offV;
            VectorXr vb = (nV>0) ? v.segment(offV, nV) : VectorXr(0);
            if (vb.size() >= 3 && reg.any_of<cardillo::C_LinearVelocity3>(e)) {
                const_cast<cardillo::C_LinearVelocity3&>(reg.get<cardillo::C_LinearVelocity3>(e)).value = vb.head<3>();
            }
            if (vb.size() >= 6 && reg.any_of<cardillo::C_AngularVelocity3>(e)) {
                const_cast<cardillo::C_AngularVelocity3&>(reg.get<cardillo::C_AngularVelocity3>(e)).value = vb.tail<3>();
            }
        }
    }
    m_world.markStateDirty();
    m_world.markForcesDirty();
}

void DynamicsAssembler::writeStateToSystem(const VectorXr& q, const VectorXr& v) {
    DynamicsAssembler::writePositionToSystem(q);
    DynamicsAssembler::writeVelocityToSystem(v);
}

void DynamicsAssembler::assignDofs() {
    auto& reg = const_cast<entt::registry&>(m_world.ecs());
    // Assign consecutive body indices to dynamic entities and compute DOF sizes in one pass
    int nextBody = 0;
    m_numQ = 0; m_numV = 0;
    auto view = reg.view<cardillo::C_PhysicsObject, cardillo::C_Position3, cardillo::C_LinearVelocity3>();
    for (auto e : view) {
        entt::entity ent = static_cast<entt::entity>(e);
        reg.emplace_or_replace<cardillo::C_BodyIndex>(ent, cardillo::C_BodyIndex{nextBody});
        ++nextBody;
        m_numQ += (index_t)m_world.getPosition(ent).size();
        m_numV += (index_t)m_world.getVelocity(ent).size();
    }
}

void DynamicsAssembler::rebuildW_() {
    const int C_all = (int)m_contacts_ptr->size();
    const int Nb = m_world.numBodies();
        // m_contact_index_orig removed; contact row mapping is stored per-contact

    // Prepare sparse triplets for W
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve((size_t)C_all * 36); // allow room for tangential rows when friction is enabled

    const auto& reg = m_world.ecs();
    const bool frictionEnabled = m_world.config().friction_enable;

    int dynContactId = 0; // index into dynamic contacts (rows in W)
    for (int i = 0; i < C_all; ++i) {
        auto& c = (*m_contacts_ptr)[i];
        const bool aDyn = reg.any_of<cardillo::C_BodyIndex>(c.a);
        const bool bDyn = reg.any_of<cardillo::C_BodyIndex>(c.b);
        // Skip static-static contacts
        if (!aDyn && !bDyn) continue;

        auto emitDirForSide = [&](const Vector3r& dir_world, const Vector3r& r_body, const Vector3r& dir_body,
                                   entt::entity ent, bool dyn, real_t s, int rowId) {
            if (!dyn) return;
            const int b = reg.get<cardillo::C_BodyIndex>(ent).b;
            if (b < 0 || b >= Nb) return;
            const int col0 = m_body_vel_offsets[(size_t)b];
            if (reg.any_of<cardillo::C_RigidBodyTag>(ent)) {
                MatrixXXr w = buildWRowRigid(dir_world, r_body, dir_body, s);
                // 6-DoF block
                for (int j = 0; j < w.cols(); ++j) {
                    real_t val = w(0,j);
                    if (val != (real_t)0) trips.emplace_back(rowId, col0 + j, val);
                }
            } else {
                MatrixXXr w = buildWRowPoint(dir_world, s);
                // 3-DoF block
                for (int j = 0; j < w.cols(); ++j) {
                    real_t val = w(0,j);
                    if (val != (real_t)0) trips.emplace_back(rowId, col0 + j, val);
                }
            }
        };

        // Row 0 for this contact: normal
        const int rowN = dynContactId;
        // record base index and default size for this contact
        c.impulse_base_index = rowN;
        c.impulse_size = 1;
        emitDirForSide(c.normal, c.pointA_body, c.normalA_body, c.a, aDyn, (real_t)-1, rowN);
        emitDirForSide(c.normal, c.pointB_body, c.normalB_body, c.b, bDyn, (real_t)+1, rowN);
        ++dynContactId;

        // Optional rows: two tangential directions if friction enabled and mu > 0
        if (frictionEnabled && c.friction_mu > (real_t)0) {
            const int rowT1 = dynContactId;
            emitDirForSide(c.tangent1, c.pointA_body, c.tangent1A_body, c.a, aDyn, (real_t)-1, rowT1);
            emitDirForSide(c.tangent1, c.pointB_body, c.tangent1B_body, c.b, bDyn, (real_t)+1, rowT1);
            ++dynContactId;

            const int rowT2 = dynContactId;
            emitDirForSide(c.tangent2, c.pointA_body, c.tangent2A_body, c.a, aDyn, (real_t)-1, rowT2);
            emitDirForSide(c.tangent2, c.pointB_body, c.tangent2B_body, c.b, bDyn, (real_t)+1, rowT2);
            ++dynContactId;
            // update impulse_size to include tangential rows
            c.impulse_size = 3;
        }
    }

    // Build W as C_dyn x totalV
    const int C_dyn = dynContactId;
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_W_sparse.resize(C_dyn, totalV);
    m_W_sparse.setFromTriplets(trips.begin(), trips.end());
    m_W_sparse.makeCompressed();
}

void DynamicsAssembler::setContactLastImpulse(int global_out_index, const cardillo::Vector3r& imp) {
    if (!m_contacts_ptr) return;
    if (global_out_index < 0 || global_out_index >= (int)m_contacts_ptr->size()) return;
    (*m_contacts_ptr)[(size_t)global_out_index].last_impulse = imp;
}

// Rebuild auxiliary block matrices derived from W and current contacts/state.
void DynamicsAssembler::rebuildInteractionW_()
{
    if (m_timings) { auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::RebuildConstraintJacobians); (void)sc; }
    // Build m_Wg/m_Wgamma and diagonals from new constraint patterns first, then legacy springs
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const auto &reg = m_world.ecs();

    std::vector<Eigen::Triplet<real_t>> tripsWg;
    std::vector<Eigen::Triplet<real_t>> tripsWgamma;
    tripsWg.reserve(1024);
    tripsWgamma.reserve(1024);

    // Per-row diagonals for C and A
    std::vector<real_t> Crows; Crows.reserve(256);
    std::vector<real_t> Arows; Arows.reserve(256);
    int springRowCounter = 0;
    int damperRowCounter = 0;
    const real_t EPS_C = (real_t)1e-10;
    const real_t EPS_A = (real_t)1e-10;

    // Emit a single 1xN row into W triplets without temporaries
    auto emitRowRef = [&](std::vector<Eigen::Triplet<real_t>>& trg,
                          int rowIndex,
                          entt::entity ent,
                          const Eigen::Ref<const Eigen::Matrix<real_t, Eigen::Dynamic, 1>>& col){
        if (!reg.any_of<cardillo::C_BodyIndex>(ent)) return;
        int b = reg.get<cardillo::C_BodyIndex>(ent).b;
        if (b < 0 || b >= (int)m_body_vel_offsets.size()-1) return;
        int row0 = m_body_vel_offsets[(size_t)b];
        int nV = m_body_vel_offsets[(size_t)b+1] - row0;
        int rows = (int)col.rows();
        int nCopy = std::min(rows, nV);
        for (int j = 0; j < nCopy; ++j) {
            real_t v = col(j);
            if (v != (real_t)0) trg.emplace_back(rowIndex, row0 + j, v);
        }
    };

    // 1) New constraint patterns (support multi-row constraints)
    const auto& patterns = m_world.constraintPatterns();
    for (const auto& uptr : patterns) {
        if (!uptr) continue;
        auto crN = uptr->getConstraint();
        const int nrows = (int)crN.Crows.size();
        // Spring rows
        for (int i = 0; i < nrows; ++i) {
            const real_t Ci = crN.Crows[i];
            if (Ci < 1 / EPS_C) {
                Crows.push_back(Ci);
                const int row = springRowCounter++;
                if (i < crN.WgA.cols()) emitRowRef(tripsWg, row, crN.a, crN.WgA.col(i));
                if (i < crN.WgB.cols()) emitRowRef(tripsWg, row, crN.b, crN.WgB.col(i));
            }
        }
        // Damper rows
        const int ndamp = (int)crN.Arows.size();
        for (int i = 0; i < ndamp; ++i) {
            const real_t Ai = crN.Arows[i];
            if (Ai < 1 / EPS_A) {
                Arows.push_back(Ai);
                const int row = damperRowCounter++;
                if (i < crN.WgammaA.cols()) emitRowRef(tripsWgamma, row, crN.a, crN.WgammaA.col(i));
                if (i < crN.WgammaB.cols()) emitRowRef(tripsWgamma, row, crN.b, crN.WgammaB.col(i));
            }
        }
    }


    // Build sparse matrices from accumulated triplets
    const int nSprings = (int)Crows.size();
    const int nDampers = (int)Arows.size();
    m_Wg.resize(nSprings, totalV);
    m_Wg.setFromTriplets(tripsWg.begin(), tripsWg.end());
    m_Wg.makeCompressed();

    m_Wgamma.resize(nDampers, totalV);
    m_Wgamma.setFromTriplets(tripsWgamma.begin(), tripsWgamma.end());
    m_Wgamma.makeCompressed();

    // Store C (per-spring) and A (per-damper) diagonals
    m_Cdiag = VectorXr::Zero((index_t)nSprings);
    for (int i = 0; i < nSprings; ++i) { m_Cdiag[i] = Crows[(size_t)i]; }

    m_Adiag = VectorXr::Zero((index_t)nDampers);
    for (int i = 0; i < nDampers; ++i) { m_Adiag[i] = Arows[(size_t)i]; }

    // Invalidate previous sparse S factorization
    m_S_sparse_lu.reset();
}

bool DynamicsAssembler::buildAndFactorS(real_t dt, real_t theta, bool includeGyroInMatrix, bool lambdaTheta)
{
    if (m_timings) { auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::BuildAndFactorS); (void)sc; }
    // Ensure current blocks are built
    rebuildInteractionW_();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const int nSprings = (int)m_Wg.rows();
    const int nDampers = (int)m_Wgamma.rows();
    const int extV = totalV + nSprings + nDampers;
    
    if (extV == 0) {
        m_S_sparse_lu.reset();
        m_S_sparse.resize(0,0);
        return true;
    }

    // Build sparse block matrix using triplets
    std::vector<Eigen::Triplet<real_t>> trips;
    const std::size_t tripEstimate = static_cast<std::size_t>(totalV)              
                                   + static_cast<std::size_t>(m_Wg.nonZeros() * 2) 
                                   + static_cast<std::size_t>(m_Wgamma.nonZeros() * 2) 
                                   + static_cast<std::size_t>(nSprings + nDampers);
    trips.reserve(tripEstimate);

    // Top-left: M diagonal
    for (int i = 0; i < totalV; ++i) {
        real_t mval = m_M_diag[i];
        if (mval != (real_t)0) trips.emplace_back(i, i, mval);
    }

    // Optionally include linearized gyroscopic term in the system matrix.
    // This adds -dt * G(omega_n) to each rigid body's rotational 3x3 block and breaks symmetry.
    if (includeGyroInMatrix) {
        const auto& reg = m_world.ecs();
        auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b < 0 || b + 1 >= (int)m_body_vel_offsets.size()) continue;
            const int off = m_body_vel_offsets[(size_t)b];
            const int nV = m_body_vel_offsets[(size_t)b + 1] - off;
            if (nV < 6) continue;
            if (!reg.all_of<cardillo::C_RigidBodyTag, cardillo::C_AngularVelocity3>(e)) continue;

            const Vector3r omega = reg.get<cardillo::C_AngularVelocity3>(e).value; // body-frame
            const Vector3r I = m_world.getInertiaDiag(e);
            const Vector3r Iomega = I.cwiseProduct(omega);
            const Matrix33r Idiag = I.asDiagonal().toDenseMatrix();
            const Matrix33r omegaSkew = skew_from_vector(omega);
            const Matrix33r IomegaSkew = skew_from_vector(Iomega);

            // G_rot = 0.5 * ( [I*omega]_x - [omega]_x * I )
            const Matrix33r Grot = (real_t)0.5 * (IomegaSkew - omegaSkew * Idiag);
            const Matrix33r corr = -dt * Grot; // M_eff = M - dt*G

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    const real_t val = corr(r, c);
                    if (val != (real_t)0) trips.emplace_back(off + 3 + r, off + 3 + c, val);
                }
            }
        }
    }

    const real_t topScale = lambdaTheta ? theta : (real_t)1.0;

    // Wg and Wgamma contributions: m_Wg is (nSprings x totalV)
    for (int k = 0; k < m_Wg.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wg, k); it; ++it) {
            int row = it.row(); // spring index
            int col = it.col(); // velocity index
            real_t v = it.value();
            // top-right: (col, totalV + row)
            trips.emplace_back(col, totalV + row, topScale * v);
            // middle-left: (totalV + row, col)
            trips.emplace_back(totalV + row, col, v);
        }
    }
    for (int k = 0; k < m_Wgamma.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wgamma, k); it; ++it) {
            int row = it.row(); int col = it.col(); real_t v = it.value();
            // top-right gamma: (col, totalV + nSprings + row)
            trips.emplace_back(col, totalV + nSprings + row, topScale * v);
            // lower-left gamma: (totalV + nSprings + row, col)
            trips.emplace_back(totalV + nSprings + row, col, v);
        }
    }

    // Middle and lower diagonal blocks (C and A terms)
    // C block over nSprings rows (assemble compliance; clamp zeros to avoid singular KKT when using SPD factorization)
    for (int i = 0; i < nSprings; ++i) {
        real_t Ci = m_Cdiag[i];
        if (!std::isfinite(Ci)) continue; // skip non-finite
        real_t cval = - (real_t)1.0 / (theta * dt * dt) * Ci;
        trips.emplace_back(totalV + i, totalV + i, cval);
    }
    // A block over nDampers rows (assemble damping compliance; clamp zeros similarly)
    for (int i = 0; i < nDampers; ++i) {
        real_t Ai = m_Adiag[i];
        if (!std::isfinite(Ai)) continue;
        real_t aval = - (real_t)1.0 / (theta * dt) * Ai;
        trips.emplace_back(totalV + nSprings + i, totalV + nSprings + i, aval);
    }

    // Build sparse matrix
    m_S_sparse.resize(extV, extV);
    m_S_sparse.setFromTriplets(trips.begin(), trips.end());
    m_S_sparse.makeCompressed();

    // Factorize using SparseLU (provide symmetry hint when applicable).
    try {
        m_S_sparse_lu.emplace();
        m_S_sparse_lu->isSymmetric(!includeGyroInMatrix && !lambdaTheta);
        m_S_sparse_lu->analyzePattern(m_S_sparse);
        m_S_sparse_lu->factorize(m_S_sparse);
        if (m_S_sparse_lu->info() != Eigen::Success) {
            m_S_sparse_lu.reset();
            std::cout << "DynamicsAssembler::buildAndFactorS: SparseLU factorization failed\n";
            return false;
        } else if (m_world.config().debug_rb) {
            std::cout << "[DynamicsAssembler] SparseLU factorization success.\n";
        }
    } catch (const std::exception& ex) {
        if (m_world.config().debug_rb) {
            std::cout << "[DynamicsAssembler] Exception during SparseLU: " << ex.what() << '\n';
        }
        m_S_sparse_lu.reset();
        return false;
    }
    return true;
}

// Build and factor S using effective mass with gyroscopic term for Stormer-Verlet
bool DynamicsAssembler::buildAndFactorS_StormerVerlet(real_t dt)
{
    if (m_timings) { auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::BuildAndFactorS); (void)sc; }
    rebuildInteractionW_();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const int nSprings = (int)m_Wg.rows();
    const int nDampers = (int)m_Wgamma.rows();
    const int extV = totalV + nSprings + nDampers;

    if (extV == 0) {
        m_S_sparse_lu.reset();
        m_S_sparse.resize(0,0);
        return true;
    }

    std::vector<Eigen::Triplet<real_t>> trips;
    const std::size_t tripEstimate = static_cast<std::size_t>(totalV)
                                   + static_cast<std::size_t>(m_Wg.nonZeros() * 2)
                                   + static_cast<std::size_t>(m_Wgamma.nonZeros() * 2)
                                   + static_cast<std::size_t>(nSprings + nDampers);
    trips.reserve(tripEstimate);

    // // Top-left: M - dt/2 * G(u) where G(u) contains gyroscopic skew term (-w_skew * I)
    // const auto &reg = m_sys.ecs();
    // auto view = reg.view<World::C_BodyIndex, World::C_PhysicsObject>();
    // for (auto [e, bi] : view.each()) {
    //     const int b = bi.b;
    //     if (b < 0 || b + 1 >= (int)m_body_vel_offsets.size()) continue;
    //     const int off = m_body_vel_offsets[(size_t)b];
    //     const int nV = m_body_vel_offsets[(size_t)b + 1] - off;

    //     // Point mass: translational block only
    //     if (reg.all_of<World::C_PointMassTag, World::C_Mass>(e)) {
    //         const real_t m = reg.get<World::C_Mass>(e).m;
    //         for (int i = 0; i < nV; ++i) trips.emplace_back(off + i, off + i, m);
    //         continue;
    //     }

    //     // Rigid body: translational mass and rotational inertia plus gyro skew term
    //     if (reg.all_of<World::C_RigidBodyTag, World::C_Mass, World::C_InertiaDiag, World::C_AngularVelocity3>(e)) {
    //         const real_t m = reg.get<World::C_Mass>(e).m;
    //         for (int i = 0; i < 3 && i < nV; ++i) trips.emplace_back(off + i, off + i, m);

    //         if (nV >= 6) {
    //             const Vector3r w = reg.get<World::C_AngularVelocity3>(e).value; // body-frame
    //             const Vector3r I = m_sys.getInertiaDiag(e);
    //             Matrix33r Idiag = I.asDiagonal().toDenseMatrix();
    //             const Vector3r Iw = I.cwiseProduct(w);
    //             auto omegaSkew = skew_from_vector(w);
    //             Matrix33r IwSkew = skew_from_vector(Iw);
    //             auto rotBlock = Idiag - ((dt * (real_t)0.5) * (IwSkew - omegaSkew * Idiag));
    //             for (int r = 0; r < 3; ++r) {
    //                 for (int c = 0; c < 3; ++c) {
    //                     const real_t val = rotBlock(r, c);
    //                     if (val != (real_t)0) trips.emplace_back(off + 3 + r, off + 3 + c, val);
    //                 }
    //             }
    //         }
    //     }
    // }

    // Top-left: M diagonal (no gyroscopic term)
    for (int i = 0; i < totalV; ++i) {
        real_t mval = m_M_diag[i];
        if (mval != (real_t)0) trips.emplace_back(i, i, mval);
    }

    // Wg and Wgamma contributions
    for (int k = 0; k < m_Wg.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wg, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            real_t v = it.value();
            trips.emplace_back(col, totalV + row, v);
            trips.emplace_back(totalV + row, col, v);
        }
    }
    for (int k = 0; k < m_Wgamma.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wgamma, k); it; ++it) {
            int row = it.row(); int col = it.col(); real_t v = it.value();
            trips.emplace_back(col, totalV + nSprings + row, v);
            trips.emplace_back(totalV + nSprings + row, col, v);
        }
    }

    // C and A blocks with Stormer-Verlet scalings
    for (int i = 0; i < nSprings; ++i) {
        real_t Ci = m_Cdiag[i];
        if (!std::isfinite(Ci)) continue;
        real_t cval = - (real_t)4.0 / (dt * dt) * Ci;
        trips.emplace_back(totalV + i, totalV + i, cval);
    }
    for (int i = 0; i < nDampers; ++i) {
        real_t Ai = m_Adiag[i];
        if (!std::isfinite(Ai)) continue;
        real_t aval = - (real_t)2.0 / dt * Ai;
        trips.emplace_back(totalV + nSprings + i, totalV + nSprings + i, aval);
    }

    m_S_sparse.resize(extV, extV);
    m_S_sparse.setFromTriplets(trips.begin(), trips.end());
    m_S_sparse.makeCompressed();

    try {
        m_S_sparse_lu.emplace();
        m_S_sparse_lu->isSymmetric(false); // gyro term breaks symmetry
        m_S_sparse_lu->analyzePattern(m_S_sparse);
        m_S_sparse_lu->factorize(m_S_sparse);
        if (m_S_sparse_lu->info() != Eigen::Success) {
            m_S_sparse_lu.reset();
            std::cout << "[DynamicsAssembler] SparseLU factorization failed\n";
            return false;
        } else if (m_world.config().debug_rb) {
            std::cout << "[DynamicsAssembler] SparseLU factorization success (Stormer-Verlet).\n";
        }
    } catch (const std::exception& ex) {
        if (m_world.config().debug_rb) {
            std::cout << "[DynamicsAssembler] Exception during SparseLU (Stormer-Verlet): " << ex.what() << '\n';
        }
        m_S_sparse_lu.reset();
        return false;
    }

    return true;
}


// Solve full extended system and return complete solution
VectorXr DynamicsAssembler::solveS(const VectorXr& rhs_ext) const
{  
    if (rhs_ext.size() == 0) return rhs_ext; // empty system: identity solve
    if (!m_S_sparse_lu.has_value()) throw std::runtime_error("DynamicsAssembler::solveS called but S matrix is not factorized");

    Eigen::VectorXd sol = m_S_sparse_lu->solve(rhs_ext);
    if (m_S_sparse_lu->info() != Eigen::Success) {
        throw std::runtime_error("DynamicsAssembler::solveS: SparseLU solve failed");
    }
    return sol;
}

void DynamicsAssembler::refreshState() {
    auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::DynamicsAssembler_RefreshState);
    bool structureChanged = false;
    if (m_world.consumeStructureDirty()) {
        structureChanged = true;
        assignDofs();
        // Recompute body offsets from ECS sizes
        const int Nb = m_world.numBodies();
        m_body_vel_offsets.assign((size_t)Nb + 1, 0);
        m_body_pos_offsets.assign((size_t)Nb + 1, 0);
        
        // First gather sizes per body index
        std::vector<int> vSizes((size_t)Nb, 0), qSizes((size_t)Nb, 0);
        const auto& reg = m_world.ecs();
        auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b; if (b < 0 || b >= Nb) continue;
            vSizes[(size_t)b] = (int)m_world.getVelocity(e).size();
            qSizes[(size_t)b] = (int)m_world.getPosition(e).size();
        }
        // Then compute prefix sums in ascending body index order
        int offV = 0, offQ = 0;
        for (int b = 0; b < Nb; ++b) {
            m_body_vel_offsets[(size_t)b] = offV;
            m_body_pos_offsets[(size_t)b] = offQ;
            offV += vSizes[(size_t)b];
            offQ += qSizes[(size_t)b];
        }
        m_body_vel_offsets[(size_t)Nb] = offV;
        m_body_pos_offsets[(size_t)Nb] = offQ;
        m_numV = offV; m_numQ = offQ;
        rebuildMass_();
        m_collision_mgr->rebuild();
    }

    if (m_world.consumeStateDirty() || structureChanged) {
        loadStateFromSystem();
    }

    if (m_world.consumeForcesDirty() || structureChanged) {
        rebuildForces_();
    }
}

void DynamicsAssembler::refreshCollisionsAndSprings(real_t dt, real_t theta, bool includeGyroInMatrix, bool lambdaTheta) {
    DerivedEntitySync::updateEntities(m_world);
    updateContactsFromSystem();
    rebuildW_();
    if (!buildAndFactorS(dt, theta, includeGyroInMatrix, lambdaTheta)) throw std::runtime_error("Failed to build and factor S matrix in DynamicsAssembler");
}

void DynamicsAssembler::refreshCollisionsAndSpringsStormerVerlet(real_t dt) {
    updateContactsFromSystem();
    rebuildW_();
    if (!buildAndFactorS_StormerVerlet(dt)) throw std::runtime_error("Failed to build and factor S matrix (Stormer-Verlet) in DynamicsAssembler");
}


} // namespace cardillo::physics