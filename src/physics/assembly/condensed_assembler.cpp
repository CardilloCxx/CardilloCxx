#include "condensed_assembler.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <utility>

#include "../../rigid_body/rigid_body.hpp"
#include "../constraints/constraints.hpp"
#include "contact_jacobian.hpp"

namespace cardillo::physics::assembly {

namespace {

// Selects the columns of `A` where `mask[i]` is true. Local copy of the identically-behaved
// helper in pgs_assembler.cpp (kept file-local to avoid a duplicate-symbol clash -- that one has
// external linkage and isn't declared in a header).
MatrixXXr selectActiveCols(const MatrixXXr& A, const std::vector<bool>& mask) {
    const int cols = (int)A.cols();
    int newCols = 0;
    for (bool b : mask)
        if (b) ++newCols;

    MatrixXXr B(A.rows(), newCols);
    int j = 0;
    for (int i = 0; i < cols; ++i) {
        if (mask[i]) B.col(j++) = A.col(i);
    }
    return B;
}

// Resolves body index + velocity-vector [offset, dof) for one side of a constraint/contact.
// Returns bodyIndex=-1, off=0, dof=0 for a static/absent entity (matches RigidBody::isStatic).
void resolveBodySide(const entt::registry& reg, const std::vector<int>& velOffsets, entt::entity e, int& bodyIndex, int& off, int& dof) {
    if (cardillo::RigidBody::isStatic(reg, e)) {
        bodyIndex = -1;
        off = 0;
        dof = 0;
        return;
    }
    bodyIndex = reg.get<cardillo::C_BodyIndex>(e).b;
    off = velOffsets[(size_t)bodyIndex];
    dof = velOffsets[(size_t)bodyIndex + 1] - off;
}

}  // namespace

CondensedTopology CondensedAssembler::buildTopology() const {
    CondensedTopology topo;
    const auto& reg = m_dyn->system().ecs();
    const auto& velOffsets = m_dyn->bodyVelOffsets();
    const auto& MinvDiag = m_dyn->MinvDiag();
    const int numBodies = velOffsets.empty() ? 0 : (int)velOffsets.size() - 1;

    std::vector<RowBlock> springBlocks, damperBlocks, contactBlocks;
    int offset = 0;

    // Springs: one pass over all constraints' active spring rows, matching PgsAssembler::Dinv()'s
    // final block order (all springs first, in constraintResults() iteration order).
    for (const auto& constraint : m_dyn->constraintResults()) {
        const auto& c_used = constraint.c_used;
        const int nSp = (int)std::count(c_used.begin(), c_used.end(), true);
        if (nSp == 0) continue;

        const bool addA = !cardillo::RigidBody::isStatic(reg, constraint.a);
        const bool addB = !cardillo::RigidBody::isStatic(reg, constraint.b);
        if (!addA && !addB) continue;

        RowBlock blk;
        blk.kind = RowBlock::Kind::Spring;
        blk.dim = nSp;
        blk.offset = offset;
        resolveBodySide(reg, velOffsets, constraint.a, blk.bodyIndexA, blk.aOff, blk.aDof);
        resolveBodySide(reg, velOffsets, constraint.b, blk.bodyIndexB, blk.bOff, blk.bDof);

        if (blk.aDof > 0) blk.Ja = selectActiveCols(constraint.WgA, c_used).transpose();  // dim x aDof
        if (blk.bDof > 0) blk.Jb = selectActiveCols(constraint.WgB, c_used).transpose();

        blk.Gii = MatrixXXr::Zero(nSp, nSp);
        if (blk.aDof > 0) blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        if (blk.bDof > 0) blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        blk.complianceDiag = VectorXr::Zero(nSp);

        springBlocks.push_back(std::move(blk));
        offset += nSp;
    }
    topo.springRows = offset;

    // Dampers: same pattern, second full pass.
    for (const auto& constraint : m_dyn->constraintResults()) {
        const auto& a_used = constraint.a_used;
        const int nDa = (int)std::count(a_used.begin(), a_used.end(), true);
        if (nDa == 0) continue;

        const bool addA = !cardillo::RigidBody::isStatic(reg, constraint.a);
        const bool addB = !cardillo::RigidBody::isStatic(reg, constraint.b);
        if (!addA && !addB) continue;

        RowBlock blk;
        blk.kind = RowBlock::Kind::Damper;
        blk.dim = nDa;
        blk.offset = offset;
        resolveBodySide(reg, velOffsets, constraint.a, blk.bodyIndexA, blk.aOff, blk.aDof);
        resolveBodySide(reg, velOffsets, constraint.b, blk.bodyIndexB, blk.bOff, blk.bDof);

        if (blk.aDof > 0) blk.Ja = selectActiveCols(constraint.WgammaA, a_used).transpose();
        if (blk.bDof > 0) blk.Jb = selectActiveCols(constraint.WgammaB, a_used).transpose();

        blk.Gii = MatrixXXr::Zero(nDa, nDa);
        if (blk.aDof > 0) blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        if (blk.bDof > 0) blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        blk.complianceDiag = VectorXr::Zero(nDa);

        damperBlocks.push_back(std::move(blk));
        offset += nDa;
    }
    topo.damperRows = offset - topo.springRows;

    // Contacts: reuse DynamicsAssembler::rebuildW_()'s own row ordering (impulse_base_index),
    // frictionless-then-frictional, instead of recomputing it.
    const auto& contacts = m_dyn->contacts();
    for (int ci = 0; ci < (int)contacts.size(); ++ci) {
        const auto& c = contacts[ci];
        if (c.impulse_base_index < 0) continue;  // static-static, skipped by rebuildW_
        const int dim = c.impulse_size;          // 1 or 3

        RowBlock blk;
        blk.kind = (dim == 1) ? RowBlock::Kind::ContactFrictionless : RowBlock::Kind::ContactFrictional;
        blk.dim = dim;
        blk.offset = topo.springRows + topo.damperRows + c.impulse_base_index;
        blk.mu = c.friction_mu;
        blk.contactIndex = ci;

        const bool aDyn = !cardillo::RigidBody::isStatic(reg, c.a);
        const bool bDyn = !cardillo::RigidBody::isStatic(reg, c.b);
        resolveBodySide(reg, velOffsets, c.a, blk.bodyIndexA, blk.aOff, blk.aDof);
        resolveBodySide(reg, velOffsets, c.b, blk.bodyIndexB, blk.bOff, blk.bDof);

        // Row order matches DynamicsAssembler::rebuildW_(): normal, then tangent1, tangent2.
        // Side A uses s=-1, side B uses s=+1 (same convention as buildContactRowByDof callers there).
        if (aDyn) {
            blk.Ja = MatrixXXr::Zero(dim, blk.aDof);
            blk.Ja.row(0) = buildContactRowByDof(blk.aDof, c.normal, c.pointA_body, c.normalA_body, (real_t)-1).transpose();
            if (dim == 3) {
                blk.Ja.row(1) = buildContactRowByDof(blk.aDof, c.tangent1, c.pointA_body, c.tangent1A_body, (real_t)-1).transpose();
                blk.Ja.row(2) = buildContactRowByDof(blk.aDof, c.tangent2, c.pointA_body, c.tangent2A_body, (real_t)-1).transpose();
            }
        }
        if (bDyn) {
            blk.Jb = MatrixXXr::Zero(dim, blk.bDof);
            blk.Jb.row(0) = buildContactRowByDof(blk.bDof, c.normal, c.pointB_body, c.normalB_body, (real_t)+1).transpose();
            if (dim == 3) {
                blk.Jb.row(1) = buildContactRowByDof(blk.bDof, c.tangent1, c.pointB_body, c.tangent1B_body, (real_t)+1).transpose();
                blk.Jb.row(2) = buildContactRowByDof(blk.bDof, c.tangent2, c.pointB_body, c.tangent2B_body, (real_t)+1).transpose();
            }
        }

        blk.Gii = MatrixXXr::Zero(dim, dim);
        if (aDyn) blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        if (bDyn) blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        blk.complianceDiag = VectorXr::Zero(dim);  // contact rows carry no compliance

        if (dim == 1)
            ++topo.frictionlessRows;
        else
            topo.frictionalRows += 3;

        contactBlocks.push_back(std::move(blk));
    }

    topo.numBilateralBlocks = (int)(springBlocks.size() + damperBlocks.size());
    topo.blocks.reserve(springBlocks.size() + damperBlocks.size() + contactBlocks.size());
    for (auto& b : springBlocks) topo.blocks.push_back(std::move(b));
    for (auto& b : damperBlocks) topo.blocks.push_back(std::move(b));
    for (auto& b : contactBlocks) topo.blocks.push_back(std::move(b));
    topo.numLambda = topo.springRows + topo.damperRows + topo.frictionlessRows + topo.frictionalRows;

    topo.blocksOfBody.assign(std::max(numBodies, 0), {});
    for (int i = 0; i < (int)topo.blocks.size(); ++i) {
        const auto& blk = topo.blocks[i];
        if (blk.bodyIndexA >= 0) topo.blocksOfBody[(size_t)blk.bodyIndexA].push_back(i);
        if (blk.bodyIndexB >= 0) topo.blocksOfBody[(size_t)blk.bodyIndexB].push_back(i);
    }

    return topo;
}

void CondensedAssembler::updateCompliance(CondensedTopology& topo, real_t dt, real_t theta) const {
    const auto& Cdiag = m_dyn->Cdiag();
    const auto& Adiag = m_dyn->Adiag();

    for (auto& blk : topo.blocks) {
        if (blk.kind == RowBlock::Kind::Spring) {
            blk.complianceDiag = Cdiag.segment(blk.offset, blk.dim) * ((real_t)1 / (theta * dt * dt));
        } else if (blk.kind == RowBlock::Kind::Damper) {
            blk.complianceDiag = Adiag.segment(blk.offset - topo.springRows, blk.dim) * ((real_t)1 / (theta * dt));
        } else {
            blk.complianceDiag = VectorXr::Zero(blk.dim);
        }

        // Match PgsAssembler::DinvDiag() -- the function PGS's own sweep loop actually calls at
        // runtime (PgsAssembler::Dinv(), the full dense block inverse, is dead code for contacts;
        // only ConjugateGradientSolver's bilateral-only path uses it). DinvDiag() is a per-row
        // DIAGONAL-only local inverse, deliberately ignoring normal/tangential (or multi-DOF)
        // coupling within a block. Using the full block inverse here changes the fixed point the
        // Gauss-Seidel/Jacobi iteration converges to on strongly-coupled multi-contact systems --
        // confirmed by direct comparison against PGS on the domino scene (diverging total contact
        // impulse and kinetic energy despite both formally "converging" by their own residual
        // check). `Gii` itself is kept as the full matrix -- the semismooth Newton local solve
        // (local_contact_newton) needs the true local coupling, that's the point of it being
        // better than this diagonal approximation.
        const VectorXr diag = blk.Gii.diagonal() + blk.complianceDiag;
        blk.GiiInv = MatrixXXr::Zero(blk.dim, blk.dim);
        for (int i = 0; i < blk.dim; ++i) blk.GiiInv(i, i) = (diag[i] > 0) ? ((real_t)1 / diag[i]) : (real_t)0;

        // // TODO: test full inversion here
        // MatrixXXr tmp = blk.Gii;
        // tmp.diagonal() += blk.complianceDiag;
        // // TOOD: Use ldlt here if implicit gyroscopic terms are not used or we can ensure that this is symmetric!
        // blk.GiiInv = tmp.partialPivLu().inverse();
    }
}

VectorXr CondensedAssembler::ufree(real_t dt, real_t theta) const {
    const auto& vn = m_dyn->vVec();
    const auto& MinvDiag = m_dyn->MinvDiag();

    VectorXr uf = vn + MinvDiag.cwiseProduct(dt * m_dyn->fVecExternal());
    if (!m_cfg.moreau_implicit_gyroscopy) uf += MinvDiag.cwiseProduct(dt * m_dyn->fVecGyroscopic());
    return uf;
}

VectorXr CondensedAssembler::rhs(const CondensedTopology& topo, real_t dt, real_t theta, const VectorXr& u_free) const {
    const auto& vn = m_dyn->vVec();
    const auto& MinvDiag = m_dyn->MinvDiag();
    const auto& MDiag = m_dyn->MDiag();
    const int totalV = m_dyn->numV();
    const int nSprings = topo.springRows;
    const int nDampers = topo.damperRows;

    VectorXr Lambda_g = m_dyn->Lambda_g();
    if ((int)Lambda_g.size() != nSprings) Lambda_g = VectorXr::Zero(nSprings);
    VectorXr Lambda_gamma = m_dyn->Lambda_gamma();
    if ((int)Lambda_gamma.size() != nDampers) Lambda_gamma = VectorXr::Zero(nDampers);

    VectorXr rhs_vel = MDiag.cwiseProduct(vn) + dt * m_dyn->fVecExternal();
    if (!m_cfg.moreau_implicit_gyroscopy)
        rhs_vel += dt * m_dyn->fVecGyroscopic();
    else
        std::cerr << "CondensedAssembler::rhs: Warning: Gyroscopic forces are treated implicitly, not implemented in condensed assembler rhs\n";

    if (m_cfg.moreau_lambda_theta && (nSprings > 0 || nDampers > 0)) {
        VectorXr corr = VectorXr::Zero(totalV);
        for (const auto& blk : topo.blocks) {
            if (blk.kind == RowBlock::Kind::Spring) {
                const VectorXr Lg = Lambda_g.segment(blk.offset, blk.dim);
                if (blk.aDof > 0) corr.segment(blk.aOff, blk.aDof).noalias() += blk.Ja.transpose() * Lg;
                if (blk.bDof > 0) corr.segment(blk.bOff, blk.bDof).noalias() += blk.Jb.transpose() * Lg;
            } else if (blk.kind == RowBlock::Kind::Damper) {
                const VectorXr Lgam = Lambda_gamma.segment(blk.offset - nSprings, blk.dim);
                if (blk.aDof > 0) corr.segment(blk.aOff, blk.aDof).noalias() += blk.Ja.transpose() * Lgam;
                if (blk.bDof > 0) corr.segment(blk.bOff, blk.bDof).noalias() += blk.Jb.transpose() * Lgam;
            }
        }
        rhs_vel -= (1.0 - theta) * corr;
    }

    const real_t beta = m_dyn->system().config().constraint_bias_factor;
    VectorXr rhs = VectorXr::Zero(topo.numLambda);

    for (const auto& blk : topo.blocks) {
        VectorXr seg = VectorXr::Zero(blk.dim);

        VectorXr WvnA_B = VectorXr::Zero(blk.dim);
        if (blk.aDof > 0) WvnA_B.noalias() += blk.Ja * vn.segment(blk.aOff, blk.aDof);
        if (blk.bDof > 0) WvnA_B.noalias() += blk.Jb * vn.segment(blk.bOff, blk.bDof);

        VectorXr WMinvRhsVel = VectorXr::Zero(blk.dim);
        if (blk.aDof > 0) WMinvRhsVel.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(rhs_vel.segment(blk.aOff, blk.aDof));
        if (blk.bDof > 0) WMinvRhsVel.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(rhs_vel.segment(blk.bOff, blk.bDof));

        if (blk.kind == RowBlock::Kind::Spring) {
            const VectorXr Crow = m_dyn->Cdiag().segment(blk.offset, blk.dim);
            const VectorXr Lg = Lambda_g.segment(blk.offset, blk.dim);

            seg = ((real_t)1 / (theta * dt * dt)) * Crow.cwiseProduct(Lg) + ((1.0 - theta) / theta) * WvnA_B + ((real_t)1 / theta) * m_dyn->C_v_vec().segment(blk.offset, blk.dim);

            if (beta > 0) {
                const VectorXr gerr = m_dyn->g_error_vec().segment(blk.offset, blk.dim);
                seg.noalias() += (-Crow.cwiseProduct(Lg) / dt + gerr) * (beta / (dt * theta));
            }
            seg.noalias() += WMinvRhsVel;
        } else if (blk.kind == RowBlock::Kind::Damper) {
            const int dOff = blk.offset - nSprings;
            seg = ((1.0 - theta) / theta) * WvnA_B + ((real_t)1 / theta) * m_dyn->A_v_vec().segment(dOff, blk.dim);
            seg.noalias() += WMinvRhsVel;
        } else {
            const int cOff = blk.offset - nSprings - nDampers;
            VectorXr ufreeTerm = VectorXr::Zero(blk.dim);
            if (blk.aDof > 0) ufreeTerm.noalias() += blk.Ja * u_free.segment(blk.aOff, blk.aDof);
            if (blk.bDof > 0) ufreeTerm.noalias() += blk.Jb * u_free.segment(blk.bOff, blk.bDof);
            seg = ufreeTerm + m_dyn->contactVVec().segment(cOff, blk.dim);
        }

        rhs.segment(blk.offset, blk.dim) = seg;
    }

    return rhs;
}

namespace {
// Returns the block's own Jacobian slice at body `b` (must be blk.bodyIndexA or blk.bodyIndexB).
const MatrixXXr& jacobianAt(const RowBlock& blk, int b) { return (blk.bodyIndexA == b) ? blk.Ja : blk.Jb; }
int offAt(const RowBlock& blk, int b) { return (blk.bodyIndexA == b) ? blk.aOff : blk.bOff; }
int dofAt(const RowBlock& blk, int b) { return (blk.bodyIndexA == b) ? blk.aDof : blk.bDof; }
}  // namespace

cardillo::misc::BlockSparseLDLT CondensedAssembler::buildBilateralFactorization(const CondensedTopology& topo) const {
    const auto& MinvDiag = m_dyn->MinvDiag();
    const int n = topo.numBilateralBlocks;

    std::vector<int> dims(n);
    std::vector<MatrixXXr> diagBlocks(n);
    for (int i = 0; i < n; ++i) {
        const auto& blk = topo.blocks[(size_t)i];
        dims[(size_t)i] = blk.dim;
        diagBlocks[(size_t)i] = blk.Gii;
        diagBlocks[(size_t)i].diagonal() += blk.complianceDiag;
    }

    // Accumulate coupling between every pair of bilateral blocks sharing a body -- a pair can in
    // principle share two bodies at once (e.g. two independent springs between the exact same two
    // bodies), so accumulate into a map keyed by the unordered pair before handing edges to
    // BlockSparseLDLT (which expects at most one edge per pair).
    std::map<std::pair<int, int>, MatrixXXr> coupling;
    for (const auto& incident : topo.blocksOfBody) {
        // Restrict to bilateral blocks incident on this body (blocksOfBody may also list contact
        // blocks, which never participate in Sbb).
        std::vector<int> bilateral;
        for (int idx : incident) {
            if (idx < n) bilateral.push_back(idx);
        }
        for (size_t ii = 0; ii < bilateral.size(); ++ii) {
            for (size_t jj = ii + 1; jj < bilateral.size(); ++jj) {
                const int p = bilateral[ii], q = bilateral[jj];
                const auto& blkP = topo.blocks[(size_t)p];
                const auto& blkQ = topo.blocks[(size_t)q];
                // The shared body: whichever of blkP's bodies also appears on blkQ's side (both
                // p and q are in this body's incident list, so it's the current `incident` body --
                // but resolve it explicitly rather than assume, in case a pair shares two bodies).
                for (int b : {blkP.bodyIndexA, blkP.bodyIndexB}) {
                    if (b < 0) continue;
                    if (b != blkQ.bodyIndexA && b != blkQ.bodyIndexB) continue;
                    const MatrixXXr& Jp = jacobianAt(blkP, b);
                    const MatrixXXr& Jq = jacobianAt(blkQ, b);
                    const int off = offAt(blkP, b), dof = dofAt(blkP, b);
                    MatrixXXr contrib = Jp * MinvDiag.segment(off, dof).asDiagonal() * Jq.transpose();  // dim_p x dim_q
                    const auto key = std::make_pair(p, q);
                    auto it = coupling.find(key);
                    if (it == coupling.end())
                        coupling.emplace(key, contrib);
                    else
                        it->second += contrib;
                }
            }
        }
    }

    std::vector<std::array<int, 2>> edgeNodes;
    std::vector<MatrixXXr> edgeBlocks;
    edgeNodes.reserve(coupling.size());
    edgeBlocks.reserve(coupling.size());
    for (auto& kv : coupling) {
        edgeNodes.push_back({kv.first.first, kv.first.second});
        edgeBlocks.push_back(std::move(kv.second));
    }

    cardillo::misc::BlockSparseLDLT ldlt;
    ldlt.build(std::move(dims), std::move(diagBlocks), edgeNodes, edgeBlocks);
    ldlt.factor();  // internal greedy minimum-degree order -- see BlockSparseLDLT's own comment
    return ldlt;
}

}  // namespace cardillo::physics::assembly
