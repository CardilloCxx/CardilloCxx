
#include "dynamics_assembler.hpp"
#include "constraints.hpp"
#include "../collision/collision_coal.hpp"
#include <Eigen/Cholesky>
#include <cmath>


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
        // One-shot: clear external force/torque components after assembling forces
        auto &reg_mut = const_cast<entt::registry&>(m_sys.ecs());
        auto viewF = reg_mut.view<PhysicsSystem::C_ExternalForce>();
        for (auto e : viewF) { reg_mut.get<PhysicsSystem::C_ExternalForce>(e).f.setZero(); }
        auto viewT = reg_mut.view<PhysicsSystem::C_ExternalTorque>();
        for (auto e : viewT) { reg_mut.get<PhysicsSystem::C_ExternalTorque>(e).tau.setZero(); }
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
        reg.emplace_or_replace<PhysicsSystem::C_BodyIndex>(ent, PhysicsSystem::C_BodyIndex{nextBody});
        ++nextBody;
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
void DynamicsAssembler::rebuildInteractionW_()
{
    // Build m_Wg/m_Wgamma and diagonals from new constraint patterns first, then legacy springs
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const auto &reg = m_sys.ecs();

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

    auto emitRowFromRow6 = [&](std::vector<Eigen::Triplet<real_t>>& trg, int rowIndex, entt::entity ent, const cardillo::Matrixr<1,6>& row){
        if (!reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(ent)) return;
        int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(ent).b;
        if (b < 0 || b >= (int)m_body_vel_offsets.size()-1) return;
        int col0 = m_body_vel_offsets[(size_t)b];
        int nV = m_body_vel_offsets[(size_t)b+1] - col0;
        if (reg.any_of<cardillo::PhysicsSystem::C_RigidBodyTag>(ent)) {
            // 6-DoF block
            int nCopy = std::min(6, nV);
            for (int j = 0; j < nCopy; ++j) {
                real_t v = row(0,j);
                if (v != (real_t)0) trg.emplace_back(rowIndex, col0 + j, v);
            }
        } else {
            // 3-DoF translational only
            int nCopy = std::min(3, nV);
            for (int j = 0; j < nCopy; ++j) {
                real_t v = row(0,j);
                if (v != (real_t)0) trg.emplace_back(rowIndex, col0 + j, v);
            }
        }
    };
    auto emitRowFromMat = [&](std::vector<Eigen::Triplet<real_t>>& trg, int rowIndex, entt::entity ent, const MatrixXXr& w){
        if (!reg.any_of<cardillo::PhysicsSystem::C_BodyIndex>(ent)) return;
        int b = reg.get<cardillo::PhysicsSystem::C_BodyIndex>(ent).b;
        if (b < 0 || b >= (int)m_body_vel_offsets.size()-1) return;
        int col0 = m_body_vel_offsets[(size_t)b];
        int nV = m_body_vel_offsets[(size_t)b+1] - col0;
        int cols = (int)w.cols();
        int nCopy = std::min(cols, nV);
        for (int j = 0; j < nCopy; ++j) {
            real_t v = w(0,j);
            if (v != (real_t)0) trg.emplace_back(rowIndex, col0 + j, v);
        }
    };

    // 1) New constraint patterns
    const auto& patterns = m_sys.constraintPatterns();
    for (const auto& uptr : patterns) {
        if (!uptr) continue;
        auto cr = uptr->getConstraint();
        // Spring part
        if (cr.C(0,0) < 1 / EPS_C) {
            Crows.push_back(cr.C(0,0));
            int row = springRowCounter++;
            emitRowFromRow6(tripsWg, row, cr.a, cr.WgA);
            emitRowFromRow6(tripsWg, row, cr.b, cr.WgB);
        }
        // Damper part
        if (cr.A(0,0) < 1 / EPS_A) {
            Arows.push_back(cr.A(0,0));
            int row = damperRowCounter++;
            emitRowFromRow6(tripsWgamma, row, cr.a, cr.WgammaA);
            emitRowFromRow6(tripsWgamma, row, cr.b, cr.WgammaB);
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
    m_S_sparse_ldlt.reset();
}

bool DynamicsAssembler::buildAndFactorS(real_t dt)
{
    // Ensure current blocks are built
    rebuildInteractionW_();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const int nSprings = (int)m_Wg.rows();
    const int nDampers = (int)m_Wgamma.rows();
    const int extV = totalV + nSprings + nDampers;

    // Build sparse block matrix using triplets
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve((size_t)totalV + (size_t)(nSprings + nDampers) * (size_t)totalV * 2 + (size_t)(nSprings + nDampers));

    // Top-left: M diagonal
    for (int i = 0; i < totalV; ++i) {
        real_t mval = m_M_diag[i];
        if (mval != (real_t)0) trips.emplace_back(i, i, mval);
    }

    // Wg and Wgamma contributions: m_Wg is (nSprings x totalV)
    for (int k = 0; k < m_Wg.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wg, k); it; ++it) {
            int row = it.row(); // spring index
            int col = it.col(); // velocity index
            real_t v = it.value();
            // top-right: (col, totalV + row)
            trips.emplace_back(col, totalV + row, v);
            // middle-left: (totalV + row, col)
            trips.emplace_back(totalV + row, col, v);
        }
    }
    for (int k = 0; k < m_Wgamma.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<real_t>::InnerIterator it(m_Wgamma, k); it; ++it) {
            int row = it.row(); int col = it.col(); real_t v = it.value();
            // top-right gamma: (col, totalV + nSprings + row)
            trips.emplace_back(col, totalV + nSprings + row, v);
            // lower-left gamma: (totalV + nSprings + row, col)
            trips.emplace_back(totalV + nSprings + row, col, v);
        }
    }

    // Middle and lower diagonal blocks (C and A terms)
    // C block over nSprings rows
    for (int i = 0; i < nSprings; ++i) {
        const real_t tiny = (real_t)1e-18;
        real_t ci = m_Cdiag[i];
        if (std::isfinite(ci) && std::abs(ci) > tiny) {
            real_t cval = - (real_t)1.0 / (dt * dt) * ci;
            if (std::isfinite(cval) && cval != (real_t)0) trips.emplace_back(totalV + i, totalV + i, cval);
        }
    }
    // A block over nDampers rows
    for (int i = 0; i < nDampers; ++i) {
        const real_t tiny = (real_t)1e-18;
        real_t ai = m_Adiag[i];
        if (std::isfinite(ai) && std::abs(ai) > tiny) {
            real_t aval = - (real_t)1.0 / dt * ai;
            if (std::isfinite(aval) && aval != (real_t)0) trips.emplace_back(totalV + nSprings + i, totalV + nSprings + i, aval);
        }
    }

    // Build sparse matrix
    m_S_sparse.resize(extV, extV);
    m_S_sparse.setFromTriplets(trips.begin(), trips.end());
    m_S_sparse.makeCompressed();

    // Factorize using SimplicialLDLT (sparse LDL^T) for symmetric matrices.
    try {

        // Create LDLT factorization object on double-valued sparse matrix
        m_S_sparse_ldlt.emplace();
        // Cast S to double for the factorization
        Eigen::SparseMatrix<double> S_double = m_S_sparse.cast<double>();
        // Todo equilibrate
        m_S_sparse_ldlt->analyzePattern(S_double);
        m_S_sparse_ldlt->factorize(S_double);
        if (m_S_sparse_ldlt->info() != Eigen::Success) {
            m_S_sparse_ldlt.reset();
            return false;
        }
    } catch (...) {
        m_S_sparse_ldlt.reset();
        return false;
    }
    return true;
}


// Solve full extended system and return complete solution
VectorXr DynamicsAssembler::solveS(const VectorXr& rhs_ext) const
{  
    if (!m_S_sparse_ldlt.has_value()) throw std::runtime_error("DynamicsAssembler::solveS called but S matrix is not factorized");

    Eigen::VectorXd sol = m_S_sparse_ldlt->solve(rhs_ext.cast<double>());
    if (m_S_sparse_ldlt->info() != Eigen::Success) {
        throw std::runtime_error("DynamicsAssembler::solveS: Sparse LDLT solve failed");
    }
    return sol.cast<real_t>();
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
        m_sys.collisionManager().rebuild();
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
