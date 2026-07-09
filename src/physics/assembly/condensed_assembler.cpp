#include "condensed_assembler.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>

#include "../../misc/block_diagonal.hpp"  // invertSmallGeneral
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

// Mirrors PjAssembler::buildAndFactorS()'s implicit-gyroscopic correction (pj_assembler.cpp:26-53)
// exactly -- same formula (Grot = 0.5*([I*omega]_x - [omega]_x*I), M_eff_rot = I - dt*Grot), same
// per-body filter (rigid body, nV>=6, has C_AngularVelocity3) -- but produces the *inverse*
// effective mass block per qualifying body instead of folding a correction into a big sparse system
// matrix: condensed's Schur-complement architecture needs an explicit Minv wherever PJ's
// direct-solve architecture needs M. Returns an empty map whenever `implicitGyro` is false (the
// overwhelming common case) -- callers treat "body absent from this map" as "use the existing plain
// -diagonal MinvDiag", so an empty map is exactly behaviorally equivalent to this feature not
// existing at all.
std::unordered_map<int, Matrixr<6, 6>> computeGyroscopicMinvBlocks(const cardillo::physics::DynamicsAssembler& dyn, real_t dt, bool implicitGyro) {
    std::unordered_map<int, Matrixr<6, 6>> blocks;
    if (!implicitGyro) return blocks;

    const auto& velOffsets = dyn.bodyVelOffsets();
    const auto& MinvDiag = dyn.MinvDiag();
    const auto& reg = dyn.system().ecs();
    auto view = reg.view<cardillo::C_BodyIndex, cardillo::C_PhysicsObject>();
    for (auto [e, bi] : view.each()) {
        const int b = bi.b;
        if (b < 0 || b + 1 >= (int)velOffsets.size()) continue;
        const int off = velOffsets[(size_t)b];
        const int nV = velOffsets[(size_t)b + 1] - off;
        if (nV < 6) continue;
        if (!reg.all_of<cardillo::C_RigidBodyTag, cardillo::C_AngularVelocity3>(e)) continue;

        const Vector3r omega = reg.get<cardillo::C_AngularVelocity3>(e).value;
        const Vector3r I = dyn.system().getInertiaDiag(e);
        const Vector3r Iomega = I.cwiseProduct(omega);
        const Matrix33r Idiag = I.asDiagonal().toDenseMatrix();
        const Matrix33r omegaSkew = cardillo::skew_from_vector(omega);
        const Matrix33r IomegaSkew = cardillo::skew_from_vector(Iomega);

        const Matrix33r Grot = (real_t)0.5 * (IomegaSkew - omegaSkew * Idiag);
        const Matrix33r Mrot = Idiag - dt * Grot;

        Matrixr<6, 6> block = Matrixr<6, 6>::Zero();
        block.topLeftCorner(3, 3) = MinvDiag.segment(off, 3).asDiagonal().toDenseMatrix();
        block.block(3, 3, 3, 3) = invertSmallGeneral(MatrixXXr(Mrot));
        blocks.emplace(b, block);
    }
    return blocks;
}

// nullptr for every body without an active override -- the overwhelming common case, so every call
// site below keeps its original plain-diagonal-MinvDiag expression unchanged in that branch (see
// each call site's comment) rather than routing through one generic "maybe non-diagonal" helper.
const Matrixr<6, 6>* gyroBlockFor(int bodyIndex, const std::unordered_map<int, Matrixr<6, 6>>& gyroBlocks) {
    if (bodyIndex < 0 || gyroBlocks.empty()) return nullptr;
    auto it = gyroBlocks.find(bodyIndex);
    return it == gyroBlocks.end() ? nullptr : &it->second;
}

}  // namespace

CondensedTopology CondensedAssembler::buildTopology(real_t dt) const {
    CondensedTopology topo;
    const auto& reg = m_dyn->system().ecs();
    const auto& velOffsets = m_dyn->bodyVelOffsets();
    const auto& MinvDiag = m_dyn->MinvDiag();
    const int numBodies = velOffsets.empty() ? 0 : (int)velOffsets.size() - 1;

    topo.gyroMinvBlocks = computeGyroscopicMinvBlocks(*m_dyn, dt, m_cfg.moreau_implicit_gyroscopy);
    const auto& gyroBlocks = topo.gyroMinvBlocks;

    std::vector<RowBlock> springBlocks, damperBlocks, contactBlocks;
    // Reserve upfront to avoid repeated vector reallocation (each doubling copies/moves every
    // existing element) as these grow via push_back below -- cheap, safe, no behavior change.
    springBlocks.reserve(m_dyn->constraintResults().size());
    damperBlocks.reserve(m_dyn->constraintResults().size());
    contactBlocks.reserve(m_dyn->contacts().size());
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
        if (blk.aDof > 0) {
            if (const auto* gA = gyroBlockFor(blk.bodyIndexA, gyroBlocks))
                blk.Gii.noalias() += blk.Ja * gA->topLeftCorner(blk.aDof, blk.aDof) * blk.Ja.transpose();
            else
                blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        }
        if (blk.bDof > 0) {
            if (const auto* gB = gyroBlockFor(blk.bodyIndexB, gyroBlocks))
                blk.Gii.noalias() += blk.Jb * gB->topLeftCorner(blk.bDof, blk.bDof) * blk.Jb.transpose();
            else
                blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        }
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
        if (blk.aDof > 0) {
            if (const auto* gA = gyroBlockFor(blk.bodyIndexA, gyroBlocks))
                blk.Gii.noalias() += blk.Ja * gA->topLeftCorner(blk.aDof, blk.aDof) * blk.Ja.transpose();
            else
                blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        }
        if (blk.bDof > 0) {
            if (const auto* gB = gyroBlockFor(blk.bodyIndexB, gyroBlocks))
                blk.Gii.noalias() += blk.Jb * gB->topLeftCorner(blk.bDof, blk.bDof) * blk.Jb.transpose();
            else
                blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        }
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
        if (aDyn) {
            if (const auto* gA = gyroBlockFor(blk.bodyIndexA, gyroBlocks))
                blk.Gii.noalias() += blk.Ja * gA->topLeftCorner(blk.aDof, blk.aDof) * blk.Ja.transpose();
            else
                blk.Gii.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).asDiagonal() * blk.Ja.transpose();
        }
        if (bDyn) {
            if (const auto* gB = gyroBlockFor(blk.bodyIndexB, gyroBlocks))
                blk.Gii.noalias() += blk.Jb * gB->topLeftCorner(blk.bDof, blk.bDof) * blk.Jb.transpose();
            else
                blk.Gii.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).asDiagonal() * blk.Jb.transpose();
        }
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

VectorXr CondensedAssembler::ufree(real_t dt, real_t theta, const CondensedTopology& topo) const {
    const auto& vn = m_dyn->vVec();
    const auto& MinvDiag = m_dyn->MinvDiag();

    VectorXr uf = vn + MinvDiag.cwiseProduct(dt * m_dyn->fVecExternal());
    if (!m_cfg.moreau_implicit_gyroscopy) uf += MinvDiag.cwiseProduct(dt * m_dyn->fVecGyroscopic());

    // Implicit-gyroscopic correction: overwrite the plain-diagonal-Minv contribution just computed
    // above, for bodies with an active override. NOTE this is NOT simply "vn + Minv_eff*dt*f_ext"
    // -- that shortcut is only valid when Minv*M=I (the plain-diagonal case, where it degenerates to
    // vn). Once Minv_eff != M^{-1}, `vn` itself must be passed through Minv_eff*M first: the system
    // being solved is M_eff*v_new = M*vn + dt*f_ext + W^T*lambda (PjAssembler's own construction --
    // M_eff only appears on the LHS, `rhs`'s top block is plain M*vn+dt*f_ext, see
    // PjAssembler::rhs()), so with no constraint force at all, v_new = Minv_eff*(M*vn + dt*f_ext).
    // No explicit gyroscopic force term is added either way -- the correction lives entirely in
    // Minv_eff. Empty loop, so a byte-for-byte no-op whenever topo.gyroMinvBlocks is empty (the
    // overwhelming common case).
    const auto& velOffsets = m_dyn->bodyVelOffsets();
    const auto& MDiag = m_dyn->MDiag();
    const auto& fExt = m_dyn->fVecExternal();
    for (const auto& [b, block] : topo.gyroMinvBlocks) {
        const int off = velOffsets[(size_t)b];
        const int dof = velOffsets[(size_t)b + 1] - off;
        const VectorXr MvnPlusFext = MDiag.segment(off, dof).cwiseProduct(vn.segment(off, dof)) + dt * fExt.segment(off, dof);
        uf.segment(off, dof) = block.topLeftCorner(dof, dof) * MvnPlusFext;
    }
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

    // Implicit gyroscopic forces (moreau_implicit_gyroscopy=true) are represented entirely through
    // the effective Minv used below (WMinvRhsVel) and upstream in buildTopology()/ufree() -- not as
    // an extra term here. Matches PjAssembler::rhs() line for line: it too uses plain M_diag*vn
    // regardless of implicitGyro, adding the explicit gyroscopic force only when NOT implicit.
    VectorXr rhs_vel = MDiag.cwiseProduct(vn) + dt * m_dyn->fVecExternal();
    if (!m_cfg.moreau_implicit_gyroscopy) rhs_vel += dt * m_dyn->fVecGyroscopic();

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
        if (blk.aDof > 0) {
            if (const auto* gA = gyroBlockFor(blk.bodyIndexA, topo.gyroMinvBlocks))
                WMinvRhsVel.noalias() += blk.Ja * (gA->topLeftCorner(blk.aDof, blk.aDof) * rhs_vel.segment(blk.aOff, blk.aDof));
            else
                WMinvRhsVel.noalias() += blk.Ja * MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(rhs_vel.segment(blk.aOff, blk.aDof));
        }
        if (blk.bDof > 0) {
            if (const auto* gB = gyroBlockFor(blk.bodyIndexB, topo.gyroMinvBlocks))
                WMinvRhsVel.noalias() += blk.Jb * (gB->topLeftCorner(blk.bDof, blk.bDof) * rhs_vel.segment(blk.bOff, blk.bDof));
            else
                WMinvRhsVel.noalias() += blk.Jb * MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(rhs_vel.segment(blk.bOff, blk.bDof));
        }

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
    // BlockSparseLDLT. `couplingFwd[{p,q}]` (p<q, i.e. the pair as first encountered) is always the
    // (p,q)-direction contribution; `couplingBwd[{p,q}]` is populated ONLY when the shared body has
    // an implicit-gyroscopic override, in which case (p,q) and (q,p) are no longer transposes of
    // each other and both must be computed independently (see BlockSparseLDLT's symmetric=false
    // mode). `anyAsym` tracks whether that ever happens anywhere in this graph -- it does not
    // flicker step to step (driven by static entity composition + the moreau_implicit_gyroscopy
    // config flag, never by the current values of Gii/complianceDiag), so the resulting edge-list
    // *shape* stays cache-stable across steps exactly like the symmetric-only case always has.
    std::map<std::pair<int, int>, MatrixXXr> couplingFwd;
    std::map<std::pair<int, int>, MatrixXXr> couplingBwd;
    bool anyAsym = false;
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
                    const auto key = std::make_pair(p, q);
                    const auto* gb = gyroBlockFor(b, topo.gyroMinvBlocks);

                    if (gb) {
                        anyAsym = true;
                        const MatrixXXr Bblk = gb->topLeftCorner(dof, dof);
                        MatrixXXr contribFwd = Jp * Bblk * Jq.transpose();  // dim_p x dim_q
                        MatrixXXr contribBwd = Jq * Bblk * Jp.transpose();  // dim_q x dim_p, NOT contribFwd^T in general
                        auto itF = couplingFwd.find(key);
                        if (itF == couplingFwd.end())
                            couplingFwd.emplace(key, contribFwd);
                        else
                            itF->second += contribFwd;
                        auto itB = couplingBwd.find(key);
                        if (itB == couplingBwd.end())
                            couplingBwd.emplace(key, contribBwd);
                        else
                            itB->second += contribBwd;
                    } else {
                        MatrixXXr contrib = Jp * MinvDiag.segment(off, dof).asDiagonal() * Jq.transpose();  // dim_p x dim_q
                        auto it = couplingFwd.find(key);
                        if (it == couplingFwd.end())
                            couplingFwd.emplace(key, contrib);
                        else
                            it->second += contrib;
                    }
                }
            }
        }
    }

    std::vector<std::array<int, 2>> edgeNodes;
    std::vector<MatrixXXr> edgeBlocks;
    if (!anyAsym) {
        // Exactly the original behavior: one undirected edge per pair, symmetric=true (default) --
        // zero cost/formula change from before this feature existed.
        edgeNodes.reserve(couplingFwd.size());
        edgeBlocks.reserve(couplingFwd.size());
        for (auto& kv : couplingFwd) {
            edgeNodes.push_back({kv.first.first, kv.first.second});
            edgeBlocks.push_back(std::move(kv.second));
        }
    } else {
        // At least one bilateral pair has a genuinely non-symmetric coupling (an implicit-
        // gyroscopic body sits between them) -- BlockSparseLDLT's symmetric=false mode needs BOTH
        // directions of EVERY edge (an edge supplied in only one direction is treated as exactly
        // zero there, not "derive the other direction"), so every still-symmetric pair in this same
        // graph also needs its true transpose supplied explicitly as its (q,p) entry.
        edgeNodes.reserve(2 * couplingFwd.size());
        edgeBlocks.reserve(2 * couplingFwd.size());
        for (auto& kv : couplingFwd) {
            const int p = kv.first.first, q = kv.first.second;
            edgeNodes.push_back({p, q});
            edgeBlocks.push_back(kv.second);
            auto itB = couplingBwd.find(kv.first);
            edgeNodes.push_back({q, p});
            edgeBlocks.push_back(itB != couplingBwd.end() ? itB->second : MatrixXXr(kv.second.transpose()));
        }
    }

    cardillo::misc::BlockSparseLDLT ldlt;
    // `dims` is passed by value (copied, not moved) so it's still available below for the cache
    // comparison -- cheap, it's just n ints.
    ldlt.build(dims, std::move(diagBlocks), edgeNodes, edgeBlocks, /*symmetric=*/!anyAsym);

    // Structural cache hit: same bilateral graph (dims + which pairs of blocks are coupled) as last
    // call -- true every step for every current scene, since springs/dampers aren't created or
    // destroyed at runtime, only their Gii/complianceDiag *values* change with dt/theta/state (both
    // already baked into `diagBlocks`/`edgeBlocks` above, independent of this check). Skip the
    // O(n^2) minimum-degree symbolic pass and reuse the cached order; factorWithOrder() still does
    // the full *numeric* factorization on the current values, so this never risks staleness in the
    // actual factorization, only in which (still-correct-for-any-order) elimination order is used.
    const bool topologyMatches = m_hasCachedBilateralOrder && dims == m_cachedBilateralDims && edgeNodes == m_cachedBilateralEdgeNodes;
    // Opt-in diagnostic (same convention as CARDILLO_DUMP_STATE, see examples/example_main.hpp):
    // prints a HIT/MISS line + running counts every call, so this cache's actual behavior on a
    // scene can be checked directly instead of assumed. Confirmed on wilberforce (999 hits / 1
    // miss over 1000 steps -- misses only the first call, as expected for a scene whose bilateral
    // topology never changes) and on dynamic_constraint (misses exactly at step 1 and at the
    // step where the second spring is added, hits every other call).
    if (getenv("CARDILLO_DEBUG_SCHUR_CACHE")) {
        static int hits = 0, misses = 0;
        (topologyMatches ? hits : misses)++;
        std::cout << "[SchurCache] " << (topologyMatches ? "HIT" : "MISS") << " n=" << n << " hits=" << hits << " misses=" << misses << std::endl;
    }
    if (topologyMatches) {
        ldlt.factorWithOrder(m_cachedBilateralOrder);
    } else {
        ldlt.factor();  // internal greedy minimum-degree order -- see BlockSparseLDLT's own comment
        m_cachedBilateralDims = std::move(dims);
        m_cachedBilateralEdgeNodes = std::move(edgeNodes);
        m_cachedBilateralOrder = ldlt.order();
        m_hasCachedBilateralOrder = true;
    }
    return ldlt;
}

}  // namespace cardillo::physics::assembly
