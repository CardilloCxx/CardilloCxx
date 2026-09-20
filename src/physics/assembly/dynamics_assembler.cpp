
#include "dynamics_assembler.hpp"
#include <Eigen/Cholesky>
#include <cmath>
#include <iostream>
#include "../../collision/collision_coal.hpp"
#include "../constraints/constraints.hpp"
#include "contact_jacobian.hpp"

namespace cardillo::physics {

void DynamicsAssembler::updateContactsFromSystem() {
    if (!m_collision_mgr) throw std::runtime_error("DynamicsAssembler::updateContactsFromSystem: no CollisionCoal provided");
    if (m_world.consumeStructureDirty()) m_collision_mgr->rebuild();
    m_collision_mgr->applyTransforms();
    m_contacts_ptr = &m_collision_mgr->detectAll();
}

void DynamicsAssembler::setLambda_g(const VectorXr& lam) { m_Lambda_g = lam; }

// ---------- Cached API (block-based) ----------

const VectorXr& DynamicsAssembler::qVec() {
    return m_q_vec;
}
const VectorXr& DynamicsAssembler::vVec() {
    return m_v_vec;
}
const VectorXr& DynamicsAssembler::fVec() {
    return m_f_vec;
}
const VectorXr& DynamicsAssembler::fVecExternal() {
    return m_f_vec_external;
}
const VectorXr& DynamicsAssembler::fVecGyroscopic() {
    return m_f_vec_gyroscopic;
}

// ---------- Rebuild helpers ----------

void DynamicsAssembler::rebuildMass_() {
    const int Nb = m_world.numBodies();
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_Minv_diag = VectorXr::Zero(totalV);
    m_M_diag = VectorXr::Zero(totalV);

    const auto& reg = m_world.ecs();
    auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b < 0 || b >= Nb) continue;
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
    const auto& reg = m_world.ecs();
    auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < Nb) {
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
    auto& reg_mut = const_cast<entt::registry&>(m_world.ecs());
    auto viewF = reg_mut.view<C_ExternalForce>();
    for (auto e : viewF) {
        reg_mut.get<C_ExternalForce>(e).f.setZero();
    }
    auto viewT = reg_mut.view<C_ExternalTorque>();
    for (auto e : viewT) {
        reg_mut.get<C_ExternalTorque>(e).tau.setZero();
    }
}

void DynamicsAssembler::loadStateFromSystem() {
    const int Nb = m_world.numBodies();
    const int totalQ = (m_body_pos_offsets.empty() ? 0 : m_body_pos_offsets.back());
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    m_q_vec = VectorXr::Zero(totalQ);
    m_v_vec = VectorXr::Zero(totalV);
    const auto& reg = m_world.ecs();
    auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
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
        }
    }
}

void DynamicsAssembler::writePositionToSystem(const VectorXr& q) {
    const auto& reg = m_world.ecs();
    auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < (int)m_body_pos_offsets.size() - 1) {
            const int offQ = m_body_pos_offsets[(size_t)b];
            const int nQ = m_body_pos_offsets[(size_t)b + 1] - offQ;
            VectorXr qb = (nQ > 0) ? q.segment(offQ, nQ) : VectorXr(0);
            if (qb.size() >= 3 && reg.any_of<C_Position3>(e)) {
                const_cast<C_Position3&>(reg.get<C_Position3>(e)).value = qb.head<3>();
            }
            if (qb.size() >= 7 && reg.any_of<C_Orientation>(e)) {
                Quaternion4r qn(qb.tail<4>());
                const Quaternion4r q_ref = reg.get<C_Orientation>(e).value;
                const_cast<C_Orientation&>(reg.get<C_Orientation>(e)).setValue(MathHelper::alignQuaternionTo(qn, q_ref));
            }
        }
    }
    m_world.markStateDirty();
    m_world.markForcesDirty();
}

void DynamicsAssembler::writeVelocityToSystem(const VectorXr& v, real_t dt) {
    auto& reg = m_world.ecs();
    auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b >= 0 && b < (int)m_body_vel_offsets.size() - 1) {
            const int offV = m_body_vel_offsets[(size_t)b];
            const int nV = m_body_vel_offsets[(size_t)b + 1] - offV;
            VectorXr vb = (nV > 0) ? v.segment(offV, nV) : VectorXr(0);
            if (vb.size() >= 3 && reg.any_of<C_LinearVelocity3>(e)) {
                auto& vlinComp = reg.get<C_LinearVelocity3>(e);
                const Vector3r prev = vlinComp.value;
                const Vector3r curr = vb.head<3>();
                vlinComp.value = curr;

                if (dt > (real_t)0) {
                    const Vector3r a = (curr - prev) / dt;
                    if (reg.any_of<C_LinearAcceleration3>(e)) {
                        reg.get<C_LinearAcceleration3>(e).value = a;
                    } else {
                        reg.emplace<C_LinearAcceleration3>(e, a);
                    }
                }
            }
            if (vb.size() >= 6 && reg.any_of<C_AngularVelocity3>(e)) {
                auto& omegaComp = reg.get<C_AngularVelocity3>(e);
                const Vector3r prev = omegaComp.value;
                const Vector3r curr = vb.tail<3>();
                omegaComp.value = curr;

                if (dt > (real_t)0) {
                    const Vector3r alpha = (curr - prev) / dt;
                    if (reg.any_of<C_AngularAcceleration3>(e)) {
                        reg.get<C_AngularAcceleration3>(e).value = alpha;
                    } else {
                        reg.emplace<C_AngularAcceleration3>(e, alpha);
                    }
                }
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
    m_numQ = 0;
    m_numV = 0;
    auto view = reg.view<C_PhysicsObject, C_Position3, C_LinearVelocity3>();
    for (auto e : view) {
        entt::entity ent = static_cast<entt::entity>(e);
        reg.emplace_or_replace<C_BodyIndex>(ent, C_BodyIndex{nextBody});
        ++nextBody;
        m_numQ += (index_t)m_world.getPosition(ent).size();
        m_numV += (index_t)m_world.getVelocity(ent).size();
    }
}

void DynamicsAssembler::rebuildW_() {
    auto sc = m_timings->scope(misc::TimingManager::TimerId::RebuildContactJacobians);
 
    auto& contacts = *m_contacts_ptr;
    const int C_all = (int)contacts.size();
    const auto& reg = m_world.ecs();
    const bool frictionEnabled = m_world.config().friction_enable;
    const auto inertial = RigidBody::RigidState::inertial();
 
    const auto sideIsDynamic = [&](const collision::ContactSide& side) {
        for (int k = 0; k < side.count; ++k) {
            if (!RigidBody::isStatic(reg, side.attachments[(size_t)k].entity)) return true;
        }
        return false;
    };
 
    // Most dissipative restitution over all attachments of a side (missing component -> defaults).
    struct Restitution {
        real_t normal, tangential;
    };
    const auto sideRestitution = [&](const collision::ContactSide& side) {
        Restitution out{std::numeric_limits<real_t>::infinity(), std::numeric_limits<real_t>::infinity()};
        for (int k = 0; k < side.count; ++k) {
            const auto* r = reg.try_get<C_Restitution>(side.attachments[(size_t)k].entity);
            const real_t n = r ? r->normal : m_cfg.restitution_default_normal;
            const real_t t = r ? r->tangential : m_cfg.restitution_default_tangential;
            out.normal = std::min(out.normal, std::max<real_t>((real_t)0, n));
            out.tangential = std::min(out.tangential, std::max<real_t>((real_t)0, t));
        }
        return out;
    };
 
    
    std::size_t attachmentCount = 0;
    for (const auto& c : contacts) attachmentCount += (std::size_t)(c.sideA.count + c.sideB.count);
    std::vector<Eigen::Triplet<real_t>> trips;
    trips.reserve(attachmentCount * (frictionEnabled ? 18 : 6));
 
    // Upper bound of 3 rows per contact; shrunk to the real row count at the end.
    m_contact_v_vec = VectorXr::Zero((index_t)C_all * 3);  // velocity contribution of static attachments
    m_mu_vec = VectorXr::Zero((index_t)C_all * 3);
    m_restitution_vec = VectorXr::Zero((index_t)C_all * 3);
 
    int dynContactId = 0;  // next free row of W
    m_numFrictionalContacts = 0;
    m_numFrictionlessContacts = 0;
 
    // ---------------------------------------------------------------------------------------------
    // Row assembly. Row layout: all frictionless contacts (1 row) first, then frictional ones
    // (normal + 2 tangential rows).
    // ---------------------------------------------------------------------------------------------
    for (bool frictionPass : {false, true}) {
        if (frictionPass && !frictionEnabled) break;
 
        for (auto& c : contacts) {
            const bool frictional = frictionEnabled && c.friction_mu > (real_t)0;
            if (frictional != frictionPass) continue;
 
            if (!sideIsDynamic(c.sideA) && !sideIsDynamic(c.sideB)) {
                c.impulse_base_index = -1;  // static-static: no rows, don't leave stale indices behind
                c.impulse_size = 0;
                continue;
            }
 
            const int nRows = frictional ? 3 : 1;
            const int rowBase = dynContactId;
            dynContactId += nRows;
            c.impulse_base_index = rowBase;
            c.impulse_size = nRows;
            (frictional ? m_numFrictionalContacts : m_numFrictionlessContacts)++;
 
            const Restitution ra = sideRestitution(c.sideA);
            const Restitution rb = sideRestitution(c.sideB);
            const std::array<Vector3r, 3> dirs = {c.normal, c.tangent1, c.tangent2};
            for (int r = 0; r < nRows; ++r) {
                m_restitution_vec[rowBase + r] = r == 0 ? std::min(ra.normal, rb.normal) : std::min(ra.tangential, rb.tangential);
                if (frictional) m_mu_vec[rowBase + r] = c.friction_mu;
            }
 
            // Adds `sign * weight * J_k` of every attachment k of `side` to the rows of this contact.
            const auto accumulateSide = [&](const collision::ContactSide& side, real_t sign) {
                for (int k = 0; k < side.count; ++k) {
                    const auto& att = side.attachments[(size_t)k];
                    const bool dyn = !RigidBody::isStatic(reg, att.entity);
                    const auto state = RigidBody::getState(reg, att.entity);
                    const real_t w = sign * att.weight;
 
                    int col0 = 0, dof = 0;
                    VectorXr veStatic;
                    if (dyn) {
                        assert(reg.all_of<C_BodyIndex>(att.entity));
                        const int b = reg.get<C_BodyIndex>(att.entity).b;
                        assert(b >= 0 && (size_t)b + 1 < m_body_vel_offsets.size());
                        col0 = m_body_vel_offsets[(size_t)b];
                        dof = m_body_vel_offsets[(size_t)b + 1] - col0;
                    } else {
                        veStatic = m_world.getVelocity(att.entity);
                        dof = (int)veStatic.size();
                    }
                    assert(dof >= 0 && dof <= 6);  // buildContactRowByDof returns a Vector6r
                    if (dof == 0) continue;
 
                    for (int r = 0; r < nRows; ++r) {
                        const Vector3r dirBody = transform::direction(dirs[(size_t)r], inertial, state);
                        const Vector6r row = buildContactRowByDof(dof, dirs[(size_t)r], att.point_body, dirBody, w);
                        if (dyn) {
                            for (int j = 0; j < dof; ++j) {
                                if (row[j] != (real_t)0) trips.emplace_back(rowBase + r, col0 + j, row[j]);
                            }
                        } else {
                            m_contact_v_vec[rowBase + r] += row.head(dof).dot(veStatic);
                        }
                    }
                }
            };
            accumulateSide(c.sideA, (real_t)-1);
            accumulateSide(c.sideB, (real_t)+1);
        }
    }
 
    // ---------------------------------------------------------------------------------------------
    // Build W as C_dyn x totalV
    // ---------------------------------------------------------------------------------------------
    const int C_dyn = dynContactId;
    const int totalV = m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back();
    m_W = TripletMatrix(C_dyn, totalV, std::make_shared<std::vector<Eigen::Triplet<real_t>>>(std::move(trips)));
    m_contact_v_vec.conservativeResize((index_t)C_dyn);
    m_mu_vec.conservativeResize((index_t)C_dyn);
    m_restitution_vec.conservativeResize((index_t)C_dyn);
}


void DynamicsAssembler::setContactLastImpulse(int global_out_index, const Vector3r& imp) {
    if (!m_contacts_ptr) return;
    if (global_out_index < 0 || global_out_index >= (int)m_contacts_ptr->size()) return;
    (*m_contacts_ptr)[(size_t)global_out_index].last_impulse = imp;
}

// Rebuild auxiliary block matrices derived from W and current contacts/state.
void DynamicsAssembler::rebuildInteractionW_() {
    auto sc = m_timings->scope(misc::TimingManager::TimerId::RebuildConstraintJacobians);

    // Build m_Wg/m_Wgamma and diagonals from new constraint patterns first, then legacy springs
    const int totalV = (m_body_vel_offsets.empty() ? 0 : m_body_vel_offsets.back());
    const auto& reg = m_world.ecs();

    std::vector<Eigen::Triplet<real_t>> tripsWg;
    std::vector<Eigen::Triplet<real_t>> tripsWgamma;
    tripsWg.reserve(1024);
    tripsWgamma.reserve(1024);

    // Per-row diagonals for C and A
    std::vector<real_t> Crows;
    std::vector<real_t> Arows;
    std::vector<real_t> C_vel;
    std::vector<real_t> A_vel;
    std::vector<real_t> g_error_vec;
    std::vector<real_t> bias_factor_vec;

    m_constraintResults.clear();

    int springRowCounter = 0;
    int damperRowCounter = 0;
    const real_t EPS_C = (real_t)1e-10;
    const real_t EPS_A = (real_t)1e-10;

    // Emit a single 1xN row into W triplets without temporaries
    auto emitColRef = [&](std::vector<Eigen::Triplet<real_t>>& trg, int rowIndex, entt::entity ent, const Eigen::Ref<const VectorXr>& col) {
        if (!reg.any_of<C_BodyIndex>(ent)) return false;
        int b = reg.get<C_BodyIndex>(ent).b;
        if (b < 0 || b >= (int)m_body_vel_offsets.size() - 1) return false;
        int row0 = m_body_vel_offsets[(size_t)b];
        int nV = m_body_vel_offsets[(size_t)b + 1] - row0;
        int rows = (int)col.rows();
        int nCopy = std::min(rows, nV);
        for (int j = 0; j < nCopy; ++j) {
            real_t v = col(j);
            if (v != (real_t)0) trg.emplace_back(rowIndex, row0 + j, v);
        }
        return true;
    };

    // 1) New constraint patterns (support multi-row constraints)
    const auto& patterns = m_world.constraintPatterns();
    for (const auto& uptr : patterns) {
        if (!uptr) continue;
        m_constraintResults.push_back(uptr->getConstraint());
        auto& constraint = m_constraintResults.back();
        constraint.c_used = std::vector<bool>(constraint.Crows.size(), false);
        constraint.a_used = std::vector<bool>(constraint.Arows.size(), false);

        const auto velSrc = uptr->getSource();
        const auto posError = constraint.positionError;

        auto velSpring = velSrc;
        auto velDamper = velSrc;
        const int nrows = (int)constraint.Crows.size();

        bool addA = !RigidBody::isStatic(reg, constraint.a);
        bool addB = !RigidBody::isStatic(reg, constraint.b);

        if (!addA && !addB) continue;
        if (!addA) {
            const VectorXr vA = m_world.getVelocity(constraint.a);
            velSpring += constraint.WgA.transpose() * vA;
            velDamper += constraint.WgammaA.transpose() * vA;
        }
        if (!addB) {
            const VectorXr vB = m_world.getVelocity(constraint.b);
            velSpring += constraint.WgB.transpose() * vB;
            velDamper += constraint.WgammaB.transpose() * vB;
        }

        // Spring rows
        for (int i = 0; i < nrows; ++i) {
            const real_t Ci = constraint.Crows[i];
            if (Ci < 1 / EPS_C) {
                Crows.push_back(Ci);
                C_vel.push_back(velSpring[i]);
                g_error_vec.push_back(posError[i]);
                bias_factor_vec.push_back(constraint.biasFactor.size() == constraint.Crows.size() ? constraint.biasFactor[i] : (real_t)1);
                const int row = springRowCounter++;
                constraint.c_used[i] = true;
                if (addA && i < constraint.WgA.cols()) emitColRef(tripsWg, row, constraint.a, constraint.WgA.col(i));
                if (addB && i < constraint.WgB.cols()) emitColRef(tripsWg, row, constraint.b, constraint.WgB.col(i));
            }
        }
        // Damper rows
        const int ndamp = (int)constraint.Arows.size();
        for (int i = 0; i < ndamp; ++i) {
            const real_t Ai = constraint.Arows[i];
            if (Ai < 1 / EPS_A) {
                Arows.push_back(Ai);
                A_vel.push_back(velDamper[i]);
                const int row = damperRowCounter++;
                constraint.a_used[i] = true;
                if (addA && i < constraint.WgammaA.cols()) emitColRef(tripsWgamma, row, constraint.a, constraint.WgammaA.col(i));
                if (addB && i < constraint.WgammaB.cols()) emitColRef(tripsWgamma, row, constraint.b, constraint.WgammaB.col(i));
            }
        }
    }

    // Build sparse matrices from accumulated triplets
    const int nSprings = (int)Crows.size();
    const int nDampers = (int)Arows.size();
    m_Wg = TripletMatrix(nSprings, totalV, std::make_shared<std::vector<Eigen::Triplet<real_t>>>(std::move(tripsWg)));
    m_Wgamma = TripletMatrix(nDampers, totalV, std::make_shared<std::vector<Eigen::Triplet<real_t>>>(std::move(tripsWgamma)));

    // Store C (per-spring) and A (per-damper) diagonals
    m_Cdiag = VectorXr::Zero((index_t)nSprings);
    m_C_v_vec = VectorXr::Zero((index_t)nSprings);
    for (int i = 0; i < nSprings; ++i) {
        m_Cdiag[i] = Crows[(size_t)i];
        m_C_v_vec[i] = C_vel[(size_t)i];
    }

    m_Adiag = VectorXr::Zero((index_t)nDampers);
    m_A_v_vec = VectorXr::Zero((index_t)nDampers);
    for (int i = 0; i < nDampers; ++i) {
        m_Adiag[i] = Arows[(size_t)i];
        m_A_v_vec[i] = A_vel[(size_t)i];
    }

    m_g_error_vec = VectorXr::Zero((index_t)nSprings);
    m_bias_factor_vec = VectorXr::Ones((index_t)nSprings);
    for (int i = 0; i < nSprings; ++i) {
        m_g_error_vec[i] = g_error_vec[(size_t)i];
        m_bias_factor_vec[i] = bias_factor_vec[(size_t)i];
    }
}

void DynamicsAssembler::refreshState() {
    auto sc = m_timings->scope(misc::TimingManager::TimerId::DynamicsAssembler_RefreshState);
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
        auto view = reg.view<C_BodyIndex, C_PhysicsObject>();
        for (auto [e, bi] : view.each()) {
            const int b = bi.b;
            if (b < 0 || b >= Nb) continue;
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
        m_numV = offV;
        m_numQ = offQ;
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

void DynamicsAssembler::updateStateDependentTerms(real_t dt) {
    {
        auto sc = m_timings->scope(misc::TimingManager::TimerId::UpdateEntities);
        DerivedEntitySync::updateEntities(m_world, dt);
    }
    updateContactsFromSystem();
    rebuildW_();
    rebuildInteractionW_();
}

}  // namespace cardillo::physics