
#include "dynamics_assembler.hpp"
#include "../collision/collision_coal.hpp"
#include "force_interaction.hpp"
#include <Eigen/Cholesky>


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
    auto& mgr = const_cast<cardillo::PhysicsSystem&>(m_sys).collisionManager();
    if (m_sys.consumeStructureDirty()) mgr.rebuild();
    mgr.applyTransforms();
    m_contacts = mgr.detectAll();
}

// ---------- Cached API (block-based) ----------

const VectorXr& DynamicsAssembler::qVec() { return m_q_vec; }
const VectorXr& DynamicsAssembler::vVec() { return m_v_vec; }
const VectorXr& DynamicsAssembler::fVec() { return m_f_vec; }

// ---------- Rebuild helpers ----------

void DynamicsAssembler::rebuildMass_() {
    const int Nb = m_sys.numBodies();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_Minv_diag = VectorXr::Zero(totalV);
    m_M_diag = VectorXr::Zero(totalV);

    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b; if (b < 0 || b >= Nb) continue;
        // Get inverse mass diagonal directly
        VectorXr MinvDiag = m_sys.getMassInverseDiag(e);
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
    const int Nb = m_sys.numBodies();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_f_vec = VectorXr::Zero(totalV);
    m_v_compat.assign((size_t)Nb, VectorXr());
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b; if (b >= 0 && b < Nb) {
            const VectorXr fb = m_sys.getForce(e);
            const int off = m_body_vel_offsets[(size_t)b];
            const int n = (int)fb.size();
            if (n > 0) std::copy(fb.data(), fb.data() + n, m_f_vec.data() + off);
        }
    }
}

void DynamicsAssembler::loadStateFromSystem() {
    const int Nb = m_sys.numBodies();
    const int totalQ = (m_body_pos_offsets.empty() ? 0 : m_body_pos_offsets.back());
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_q_vec = VectorXr::Zero(totalQ);
    m_v_vec = VectorXr::Zero(totalV);
    m_v_compat.assign((size_t)Nb, VectorXr()); // deprecated compatibility only
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < Nb) {
            const VectorXr qb = m_sys.getPosition(e);
            const VectorXr vb = m_sys.getVelocity(e);
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

void DynamicsAssembler::writeStateToSystem(const VectorXr& q, const VectorXr& v) {
    const auto& reg = m_sys.ecs();
    auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < (int)m_body_pos_offsets.size()-1) {
            const int offQ = m_body_pos_offsets[(size_t)b];
            const int nQ = m_body_pos_offsets[(size_t)b+1] - offQ;
            VectorXr qb = (nQ>0) ? q.segment(offQ, nQ) : VectorXr(0);
            if (qb.size() >= 3 && reg.any_of<PhysicsSystem::C_Position3>(e)) {
                const_cast<PhysicsSystem::C_Position3&>(reg.get<PhysicsSystem::C_Position3>(e)).value = qb.head<3>();
            }
            if (qb.size() >= 7 && reg.any_of<PhysicsSystem::C_Orientation>(e)) {
                Quaternion4r qn(qb.tail<4>()); qn.normalize();
                const_cast<PhysicsSystem::C_Orientation&>(reg.get<PhysicsSystem::C_Orientation>(e)).q = qn;
            }
        }
        if (b >= 0 && b < (int)m_body_vel_offsets.size()-1) {
            const int offV = m_body_vel_offsets[(size_t)b];
            const int nV = m_body_vel_offsets[(size_t)b+1] - offV;
            VectorXr vb = (nV>0) ? v.segment(offV, nV) : VectorXr(0);
            if (vb.size() >= 3 && reg.any_of<PhysicsSystem::C_LinearVelocity3>(e)) {
                const_cast<PhysicsSystem::C_LinearVelocity3&>(reg.get<PhysicsSystem::C_LinearVelocity3>(e)).value = vb.head<3>();
            }
            if (vb.size() >= 6 && reg.any_of<PhysicsSystem::C_AngularVelocity3>(e)) {
                const_cast<PhysicsSystem::C_AngularVelocity3&>(reg.get<PhysicsSystem::C_AngularVelocity3>(e)).value = vb.tail<3>();
            }
        }
    }
    m_sys.markStateDirty();
    m_sys.markForcesDirty();
}

void DynamicsAssembler::assignDofs() {
    auto& reg = const_cast<entt::registry&>(m_sys.ecs());
    // Assign consecutive body indices to dynamic entities and compute DOF sizes in one pass
    int nextBody = 0;
    m_numQ = 0; m_numV = 0;
    auto view = reg.view<PhysicsSystem::C_PhysicsObject, PhysicsSystem::C_Position3, PhysicsSystem::C_LinearVelocity3>();
    for (auto e : view) {
        entt::entity ent = static_cast<entt::entity>(e);
        if (!reg.any_of<PhysicsSystem::C_BodyIndex>(ent)) {
            reg.emplace<PhysicsSystem::C_BodyIndex>(ent, PhysicsSystem::C_BodyIndex{nextBody});
            ++nextBody;
        }
        m_numQ += (index_t)m_sys.getPosition(ent).size();
        m_numV += (index_t)m_sys.getVelocity(ent).size();
    }
}

void DynamicsAssembler::rebuildW_() {
    const int C_all = (int)m_contacts.size();
    const int Nb = m_sys.numBodies();
    m_contact_index_orig.clear();
    m_contact_index_orig.reserve((size_t)C_all);

    // Prepare sparse triplets for W
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve((size_t)C_all * 36); // allow room for tangential rows when friction is enabled

    const auto& reg = m_sys.ecs();
    const bool frictionEnabled = m_sys.config().friction_enable;

    int dynContactId = 0; // index into dynamic contacts (rows in W)
    for (int i = 0; i < C_all; ++i) {
        const auto& c = m_contacts[i];
        const bool aDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.a);
        const bool bDyn = reg.any_of<PhysicsSystem::C_BodyIndex>(c.b);
        // Skip static-static contacts
        if (!aDyn && !bDyn) continue;

        auto emitDirForSide = [&](const Vector3r& dir_world, const Vector3r& r_body, const Vector3r& dir_body,
                                   entt::entity ent, bool dyn, real_t s, int rowId) {
            if (!dyn) return;
            const int b = reg.get<PhysicsSystem::C_BodyIndex>(ent).b;
            if (b < 0 || b >= Nb) return;
            const int col0 = m_body_vel_offsets[(size_t)b];
            if (reg.any_of<PhysicsSystem::C_RigidBodyTag>(ent)) {
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
        emitDirForSide(c.normal, c.pointA_body, c.normalA_body, c.a, aDyn, (real_t)-1, rowN);
        emitDirForSide(c.normal, c.pointB_body, c.normalB_body, c.b, bDyn, (real_t)+1, rowN);
        m_contact_index_orig.push_back(i);
        ++dynContactId;

        // Optional rows: two tangential directions if friction enabled and mu > 0
        if (frictionEnabled && c.friction_mu > (real_t)0) {
            const int rowT1 = dynContactId;
            emitDirForSide(c.tangent1, c.pointA_body, c.tangent1A_body, c.a, aDyn, (real_t)-1, rowT1);
            emitDirForSide(c.tangent1, c.pointB_body, c.tangent1B_body, c.b, bDyn, (real_t)+1, rowT1);
            m_contact_index_orig.push_back(i);
            ++dynContactId;

            const int rowT2 = dynContactId;
            emitDirForSide(c.tangent2, c.pointA_body, c.tangent2A_body, c.a, aDyn, (real_t)-1, rowT2);
            emitDirForSide(c.tangent2, c.pointB_body, c.tangent2B_body, c.b, bDyn, (real_t)+1, rowT2);
            m_contact_index_orig.push_back(i);
            ++dynContactId;
        }
    }

    // Build W as C_dyn x totalV
    const int C_dyn = dynContactId;
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_W_sparse.resize(C_dyn, totalV);
    m_W_sparse.setFromTriplets(trips.begin(), trips.end());
    m_W_sparse.makeCompressed();
}

// Rebuild auxiliary block matrices derived from W and current contacts/state.
void DynamicsAssembler::rebuidInteractionW_()
{
    // Build m_Wg from force interactions (spring constraints) and assemble global K/D
    auto &mgr = const_cast<cardillo::PhysicsSystem&>(m_sys).forceManager();
    auto &constraints = mgr.constraints();
    const int numSprings = (int)constraints.size();
    const int Cdyn = 6 * numSprings; // 6 rows per spring
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());

    // Prepare triplets for Wg (Cdyn x totalV) and K/D (Cdyn x Cdyn)
    std::vector<Eigen::Triplet<real_t>> tripsWg;
    std::vector<Eigen::Triplet<real_t>> tripsWgamma;
    std::vector<Eigen::Triplet<real_t>> tripsK;
    std::vector<Eigen::Triplet<real_t>> tripsD;
    tripsWg.reserve((size_t)numSprings * 12 * 6);
    tripsWgamma.reserve((size_t)numSprings * 12 * 6);
    tripsK.reserve((size_t)numSprings * 36);
    tripsD.reserve((size_t)numSprings * 36);

    // For each constraint, build local Wg blocks and append into global matrices
    const auto &reg = m_sys.ecs();
    m_gcat = VectorXr::Zero(Cdyn);
    m_gdotcat = VectorXr::Zero(Cdyn);
    for (int i = 0; i < numSprings; ++i) {
        auto &c = constraints[(size_t)i];
        // Ensure per-constraint cached values are up-to-date (registry should already be set by passthrough)
        c.recomputeDeformations();

        Matrix6r WgA, WgB, WgammaA, WgammaB;
        c.computeWBlocks(WgA, WgB, WgammaA, WgammaB);

        const int row0 = 6 * i;

        // Store concatenated g and gdot
        m_gcat.segment<6>(row0) = c.g;
        m_gdotcat.segment<6>(row0) = c.gdot;

        // Body A
        if (c.bodyA != entt::null && reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(c.bodyA)) {
            int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(c.bodyA).b;
            if (b >= 0 && b < (int)m_body_vel_offsets.size()-1) {
                int col0 = m_body_vel_offsets[(size_t)b];
                int nV = m_body_vel_offsets[(size_t)b+1] - col0;
                int nCopy = std::min(6, nV);
                for (int r = 0; r < 6; ++r) {
                    for (int j = 0; j < nCopy; ++j) {
                        real_t val = WgA(r,j);
                        if (val != (real_t)0) tripsWg.emplace_back(row0 + r, col0 + j, val);
                    }
                }
                // Wgamma contributions for Body A
                for (int r = 0; r < 6; ++r) {
                    for (int j = 0; j < nCopy; ++j) {
                        real_t val = WgammaA(r,j);
                        if (val != (real_t)0) tripsWgamma.emplace_back(row0 + r, col0 + j, val);
                    }
                }
            }
        }

        // Body B
        if (c.bodyB != entt::null && reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(c.bodyB)) {
            int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(c.bodyB).b;
            if (b >= 0 && b < (int)m_body_vel_offsets.size()-1) {
                int col0 = m_body_vel_offsets[(size_t)b];
                int nV = m_body_vel_offsets[(size_t)b+1] - col0;
                int nCopy = std::min(6, nV);
                for (int r = 0; r < 6; ++r) {
                    for (int j = 0; j < nCopy; ++j) {
                        real_t val = WgB(r,j);
                        if (val != (real_t)0) tripsWg.emplace_back(row0 + r, col0 + j, val);
                    }
                }
                // Wgamma contributions for Body B
                for (int r = 0; r < 6; ++r) {
                    for (int j = 0; j < nCopy; ++j) {
                        real_t val = WgammaB(r,j);
                        if (val != (real_t)0) tripsWgamma.emplace_back(row0 + r, col0 + j, val);
                    }
                }
            }
        }

        // Assemble K and D block-diagonal entries for this constraint
        for (int r = 0; r < 6; ++r) {
            for (int ccol = 0; ccol < 6; ++ccol) {
                real_t kval = c.K(r, ccol);
                real_t dval = c.D(r, ccol);
                if (kval != (real_t)0) tripsK.emplace_back(row0 + r, row0 + ccol, kval);
                if (dval != (real_t)0) tripsD.emplace_back(row0 + r, row0 + ccol, dval);
            }
        }
    }

    // Build sparse matrices
    m_Wg.resize(Cdyn, totalV);
    m_Wg.setFromTriplets(tripsWg.begin(), tripsWg.end());
    m_Wg.makeCompressed();

    m_Wgamma.resize(Cdyn, totalV);
    m_Wgamma.setFromTriplets(tripsWgamma.begin(), tripsWgamma.end());
    m_Wgamma.makeCompressed();

    m_K.resize(Cdyn, Cdyn);
    m_K.setFromTriplets(tripsK.begin(), tripsK.end());
    m_K.makeCompressed();

    m_D.resize(Cdyn, Cdyn);
    m_D.setFromTriplets(tripsD.begin(), tripsD.end());
    m_D.makeCompressed();

    // Invalidate previous S factorization
    m_S_llt.reset();
}

bool DynamicsAssembler::buildAndFactorS(real_t dt)
{
    // Ensure current blocks are built
    rebuidInteractionW_();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
   
    Eigen::SparseMatrix<real_t> Msp(totalV, totalV);
    std::vector<Eigen::Triplet<real_t>> tripsM;
    tripsM.reserve((size_t)totalV);
    for (int i = 0; i < totalV; ++i) {
        real_t mval = m_M_diag[i];
        if (mval != (real_t)0) tripsM.emplace_back(i, i, mval);
    }
    Msp.setFromTriplets(tripsM.begin(), tripsM.end());

    // Diagnostics: print sizes and sparsity of contributing matrices so we can trace zero-entries
    std::cerr << "[DynamicsAssembler] S-assembly inputs: totalV=" << totalV
              << " Msp.nnz=" << Msp.nonZeros()
              << " m_M_diag_nonzero=" << (int)(std::count_if(m_M_diag.data(), m_M_diag.data()+totalV, [](real_t v){return v!=0;}))
              << " m_Wg(" << m_Wg.rows() << "," << m_Wg.cols() << ").nnz=" << m_Wg.nonZeros()
              << " m_K(" << m_K.rows() << "," << m_K.cols() << ").nnz=" << m_K.nonZeros()
              << " m_Wgamma(" << m_Wgamma.rows() << "," << m_Wgamma.cols() << ").nnz=" << m_Wgamma.nonZeros()
              << " m_D.nnz=" << m_D.nonZeros() << std::endl;

    // Build S = M + dt^2 * Wg^T * K * Wg + dt * Wgamma^T * D * Wgamma
    Eigen::SparseMatrix<real_t> Sterm1(totalV, totalV);
    Eigen::SparseMatrix<real_t> Sterm2(totalV, totalV);

    if (m_Wg.rows() > 0 && m_K.rows() > 0) {
        // compute dt^2 * (Wg^T * (K * Wg))
        Sterm1 = (dt * dt) * (m_Wg.transpose() * (m_K * m_Wg));
        std::cerr << "[DynamicsAssembler] Sterm1 built; nnz=" << Sterm1.nonZeros() << std::endl;
    }

    if (m_Wgamma.rows() > 0 && m_D.rows() > 0) {
        // compute dt * (Wgamma^T * (D * Wgamma))
        Sterm2 = dt * (m_Wgamma.transpose() * (m_D * m_Wgamma));
        std::cerr << "[DynamicsAssembler] Sterm2 built; nnz=" << Sterm2.nonZeros() << std::endl;
    }

    m_S = Msp + Sterm1 + Sterm2;
    m_S.makeCompressed();

    // Factorize S using SimplicialLLT for self-adjoint positive definite matrices
    try {
        m_S_llt.emplace();
        m_S_llt->compute(m_S);
        if (m_S_llt->info() != Eigen::Success) {
            m_S_llt.reset();
            return false;
        }
    } catch (...) {
        m_S_llt.reset();
        return false;
    }
    return true;
}

VectorXr DynamicsAssembler::solveWithS(const VectorXr& rhs) const
{
    if (!m_S_llt.has_value()) throw std::runtime_error("S factorization not available");
    Eigen::VectorXd rhsd = rhs.cast<double>();
    Eigen::VectorXd x = m_S_llt->solve(rhsd);
    VectorXr out((index_t)x.size());
    for (int i = 0; i < (int)x.size(); ++i) out[i] = (real_t)x[i];
    return out;
}

void DynamicsAssembler::refreshState() {
    bool structureChanged = false;
    if (m_sys.consumeStructureDirty()) {
        structureChanged = true;
        assignDofs();
        // Recompute body offsets from ECS sizes
        const int Nb = m_sys.numBodies();
        m_body_vel_offsets.assign((size_t)Nb + 1, 0);
        m_body_pos_offsets.assign((size_t)Nb + 1, 0);
        
        // First gather sizes per body index
        std::vector<int> vSizes((size_t)Nb, 0), qSizes((size_t)Nb, 0);
        const auto& reg = m_sys.ecs();
        auto view = reg.view<PhysicsSystem::C_BodyIndex, PhysicsSystem::C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b; if (b < 0 || b >= Nb) continue;
            vSizes[(size_t)b] = (int)m_sys.getVelocity(e).size();
            qSizes[(size_t)b] = (int)m_sys.getPosition(e).size();
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
    }

    if (m_sys.consumeStateDirty() || structureChanged) {
        loadStateFromSystem();
    }

    if (m_sys.consumeForcesDirty() || structureChanged) {
        rebuildForces_();
    }
}

void DynamicsAssembler::refreshCollisionsAndSprings(real_t dt) {
    updateContactsFromSystem();
    rebuildW_();

    if (!buildAndFactorS(dt)) throw std::runtime_error("Failed to build and factor S matrix in DynamicsAssembler");
}


} // namespace cardillo::physics
