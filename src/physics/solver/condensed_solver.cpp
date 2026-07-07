#include "condensed_solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>

#ifdef CARDILLO_HAVE_OPENMP
#include <omp.h>
#endif

#include "../../misc/graph_coloring.hpp"
#include "local_contact_newton.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

using cardillo::physics::assembly::CondensedTopology;
using cardillo::physics::assembly::RowBlock;

namespace {

// Same clip as ProjectedGaussSeidel's project()/project_all(): normal impulse (index 0) clamped to
// <= 0, tangential (indices 1,2) scaled toward the Coulomb disk of radius -mu*lambda[0]. No-op for
// Spring/Damper blocks (equality constraints, unconstrained).
void projectBlock(RowBlock::Kind kind, real_t mu, VectorXr& lam) {
    if (kind == RowBlock::Kind::Spring || kind == RowBlock::Kind::Damper) return;

    lam[0] = std::min(lam[0], (real_t)0);
    if (kind == RowBlock::Kind::ContactFrictionless) return;

    const real_t tnorm = std::sqrt(lam[1] * lam[1] + lam[2] * lam[2]);
    if (tnorm <= (real_t)0) return;
    const real_t s = std::min((real_t)1, (-mu * lam[0]) / tnorm);
    lam[1] *= s;
    lam[2] *= s;
}

// One block's local residual r = rhs_block - Ja*u_corr_A - Jb*u_corr_B - complianceDiag*lambda_block.
VectorXr blockResidual(const RowBlock& blk, const VectorXr& rhs, const VectorXr& u_corr, const VectorXr& lambda) {
    VectorXr r = rhs.segment(blk.offset, blk.dim);
    if (blk.aDof > 0) r.noalias() -= blk.Ja * u_corr.segment(blk.aOff, blk.aDof);
    if (blk.bDof > 0) r.noalias() -= blk.Jb * u_corr.segment(blk.bOff, blk.bDof);
    r.noalias() -= blk.complianceDiag.cwiseProduct(lambda.segment(blk.offset, blk.dim));
    return r;
}

// Applies one block's impulse delta to u_corr (in place): u_corr += Minv * J^T * dlambda.
void scatterDelta(const RowBlock& blk, const VectorXr& MinvDiag, const VectorXr& dlambda, VectorXr& u_corr) {
    if (blk.aDof > 0) u_corr.segment(blk.aOff, blk.aDof).noalias() += MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(blk.Ja.transpose() * dlambda);
    if (blk.bDof > 0) u_corr.segment(blk.bOff, blk.bDof).noalias() += MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(blk.Jb.transpose() * dlambda);
}

bool isContactKind(RowBlock::Kind kind) { return kind == RowBlock::Kind::ContactFrictionless || kind == RowBlock::Kind::ContactFrictional; }

// Shared local-solve step used by every sweep driver: for frictional contact blocks when
// `useNewton` is set, tries the guarded Alart-Curnier semismooth Newton solve first (see
// local_contact_newton.hpp) -- it either converges to a locally-exact solution in a handful of
// iterations, or reports failure (singular Jacobian, non-finite step, non-decreasing residual, or
// no convergence within its iteration cap), in which case we fall through to the always-available
// linear step + projection below. Returns the new (projected/converged) impulse for this block;
// caller computes the delta.
VectorXr computeBlockUpdate(const RowBlock& blk, const VectorXr& r, real_t relaxation, real_t alpha, const VectorXr& lamOld, bool useNewton, const NewtonACParams& newtonParams) {
    if (useNewton && blk.kind == RowBlock::Kind::ContactFrictional) {
        Vector3r lam3 = lamOld;
        if (solveContactBlockNewtonAC(blk.Gii, r, blk.mu, lam3, newtonParams)) return VectorXr(lam3);
    }

    VectorXr dlam = relaxation * (blk.GiiInv * r);
    if (isContactKind(blk.kind)) dlam *= alpha;
    VectorXr lamNew = lamOld + dlam;
    projectBlock(blk.kind, blk.mu, lamNew);
    return lamNew;
}

// Sequential block-Gauss-Seidel sweep: propagates u_corr immediately after each block, matching
// ProjectedGaussSeidelSolver's own loop structure exactly, just expressed over RowBlocks.
void gaussSeidelSweep(const CondensedTopology& topo, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton, const NewtonACParams& newtonParams,
                       VectorXr& lambda, VectorXr& u_corr) {
    for (const auto& blk : topo.blocks) {
        const VectorXr r = blockResidual(blk, rhs, u_corr, lambda);
        const VectorXr lamOld = lambda.segment(blk.offset, blk.dim);
        const VectorXr lamNew = computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams);

        scatterDelta(blk, MinvDiag, lamNew - lamOld, u_corr);
        lambda.segment(blk.offset, blk.dim) = lamNew;
    }
}

// Builds the block adjacency graph for coloring: two blocks are adjacent iff they share a body
// (i.e. both appear in that body's incident-block list from CondensedTopology::blocksOfBody).
// Deliberately NOT the dense Schur-complement/Delassus sparsity pattern (which for a compliant
// chain would be fully dense and collapse coloring into one color) -- this is the same
// shared-body adjacency Vivace/Bullet-style contact coloring uses.
std::vector<std::vector<int>> buildBlockAdjacency(const CondensedTopology& topo) {
    std::vector<std::vector<int>> adj(topo.blocks.size());
    for (const auto& incident : topo.blocksOfBody) {
        for (size_t i = 0; i < incident.size(); ++i) {
            for (size_t j = i + 1; j < incident.size(); ++j) {
                adj[(size_t)incident[i]].push_back(incident[j]);
                adj[(size_t)incident[j]].push_back(incident[i]);
            }
        }
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return adj;
}

// Graph-colored Gauss-Seidel: sequential across colors, (eventually) parallel within a color.
// Safe by construction -- no two blocks in the same color class share a body, so their reads/
// writes of u_corr touch disjoint index ranges.
void coloredSweep(const CondensedTopology& topo, const cardillo::misc::Coloring& coloring, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                   const NewtonACParams& newtonParams, VectorXr& lambda, VectorXr& u_corr) {
    const int numColors = coloring.numColors;
    // Single team for the whole sweep (all colors), not one team per color: spawning/joining an
    // OpenMP team has real overhead, and a sweep can have dozens of colors -- doing that per color,
    // every sweep, every solve() call, dominates wall-clock on scenes with many contacts. Each
    // `#pragma omp for` still has its own implicit barrier, so the color-sequencing invariant
    // (no color starts until the previous one's writes are all visible) is unchanged.
#pragma omp parallel
    for (int c = 0; c < numColors; ++c) {
        const auto& cls = coloring.colorClasses[(size_t)c];
        const int n = (int)cls.size();
#pragma omp for schedule(dynamic)
        for (int k = 0; k < n; ++k) {
            const int idx = cls[(size_t)k];
            const auto& blk = topo.blocks[(size_t)idx];
            const VectorXr r = blockResidual(blk, rhs, u_corr, lambda);
            const VectorXr lamOld = lambda.segment(blk.offset, blk.dim);
            const VectorXr lamNew = computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams);

            // Safe: no two blocks in `cls` share a body (that's the coloring invariant), so these
            // writes land in disjoint slices of u_corr/lambda across threads -- no race, no atomics.
            scatterDelta(blk, MinvDiag, lamNew - lamOld, u_corr);
            lambda.segment(blk.offset, blk.dim) = lamNew;
        }
        // implicit barrier at the end of `omp for` here, before the next color starts
    }
}

// Jacobi sweep: two passes, no shared-state races even when parallelized (see design note in
// condensed_solver.hpp). Pass 1 updates each block's own impulse against a frozen snapshot of
// u_corr from the previous sweep. Pass 2 has each BODY gather corrections from only its own
// incident blocks -- never a block scattering into shared state -- so the reduction needs no
// synchronization at all, not even atomics.
void jacobiSweep(const CondensedTopology& topo, const std::vector<int>& bodyVelOffsets, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                  const NewtonACParams& newtonParams, VectorXr& lambda, VectorXr& u_corr) {
    const VectorXr u_read = u_corr;
    const int numBlocks = (int)topo.blocks.size();
    const int numBodies = (int)topo.blocksOfBody.size();
    VectorXr dlambda = VectorXr::Zero(lambda.size());
    VectorXr u_write = u_read;

    // Single team for both passes (one spawn/join per sweep, not two): each `#pragma omp for` has
    // its own implicit barrier, so pass 2 still never starts until every thread has finished
    // writing pass 1's dlambda -- the correctness argument is unchanged, only the team lifetime.
#pragma omp parallel
    {
        // Pass 1: each thread handles a disjoint set of block indices, writing only to that
        // block's own (non-overlapping) segment of lambda/dlambda -- no aliasing across threads.
#pragma omp for schedule(static)
        for (int i = 0; i < numBlocks; ++i) {
            const auto& blk = topo.blocks[i];
            const VectorXr r = blockResidual(blk, rhs, u_read, lambda);
            const VectorXr lamOld = lambda.segment(blk.offset, blk.dim);
            const VectorXr lamNew = computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams);

            dlambda.segment(blk.offset, blk.dim) = lamNew - lamOld;
            lambda.segment(blk.offset, blk.dim) = lamNew;
        }
        // implicit barrier here, before pass 2 starts

        // Pass 2: each thread handles a disjoint set of BODIES, gathering only from that body's
        // own incident blocks and writing only to that body's own segment of u_write -- again no
        // aliasing across threads, and no atomics needed (gather, not scatter).
#pragma omp for schedule(static)
        for (int b = 0; b < numBodies; ++b) {
            const auto& incident = topo.blocksOfBody[b];
            if (incident.empty()) continue;
            const int off = bodyVelOffsets[(size_t)b];
            const int dof = bodyVelOffsets[(size_t)b + 1] - off;

            VectorXr acc = VectorXr::Zero(dof);
            for (int blockIdx : incident) {
                const auto& blk = topo.blocks[blockIdx];
                const VectorXr d = dlambda.segment(blk.offset, blk.dim);
                if (blk.bodyIndexA == b) acc.noalias() += blk.Ja.transpose() * d;
                if (blk.bodyIndexB == b) acc.noalias() += blk.Jb.transpose() * d;
            }
            u_write.segment(off, dof) = u_read.segment(off, dof) + MinvDiag.segment(off, dof).cwiseProduct(acc);
        }
    }
    u_corr = std::move(u_write);
}

// Chaotic / randomized Gauss-Seidel: deliberately unsynchronized updates in (periodically
// reshuffled) random order. `order` is a permutation of block indices, so within one call each
// block index is touched by exactly one thread -- writes to `lambda` are race-free by the same
// disjoint-segment argument as every other sweep. What's genuinely different here is `u_corr`:
// unlike Jacobi's frozen snapshot or colored's provably-disjoint writes, here two blocks sharing a
// body CAN run concurrently and race on that body's accumulator -- that staleness is the point
// (Chazan-Miranker chaotic relaxation), but it must be implemented with real atomics, not a data
// race: plain-`double` unsynchronized read-modify-write is undefined behavior in C++, not just
// imprecise. `real_t` is `double`; C++20 guarantees `std::atomic<double>::fetch_add`/`load` exist
// (P0020R6), so no hand-rolled CAS loop is needed.
void chaoticSweep(const CondensedTopology& topo, std::vector<int>& order, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                   const NewtonACParams& newtonParams, int reshuffleInterval, int iter, std::mt19937& rng, VectorXr& lambda, std::vector<std::atomic<real_t>>& u_corr_atomic) {
    if (iter % reshuffleInterval == 0) std::shuffle(order.begin(), order.end(), rng);

    const int n = (int)order.size();
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < n; ++k) {
        const int idx = order[(size_t)k];
        const auto& blk = topo.blocks[(size_t)idx];

        VectorXr uA(blk.aDof), uB(blk.bDof);
        for (int d = 0; d < blk.aDof; ++d) uA[d] = u_corr_atomic[(size_t)(blk.aOff + d)].load(std::memory_order_relaxed);
        for (int d = 0; d < blk.bDof; ++d) uB[d] = u_corr_atomic[(size_t)(blk.bOff + d)].load(std::memory_order_relaxed);

        VectorXr r = rhs.segment(blk.offset, blk.dim);
        if (blk.aDof > 0) r.noalias() -= blk.Ja * uA;
        if (blk.bDof > 0) r.noalias() -= blk.Jb * uB;
        r.noalias() -= blk.complianceDiag.cwiseProduct(lambda.segment(blk.offset, blk.dim));

        const VectorXr lamOld = lambda.segment(blk.offset, blk.dim);
        const VectorXr lamNew = computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams);
        lambda.segment(blk.offset, blk.dim) = lamNew;  // safe: `order` is a permutation this sweep

        const VectorXr d = lamNew - lamOld;
        if (blk.aDof > 0) {
            const VectorXr incA = MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(blk.Ja.transpose() * d);
            for (int dd = 0; dd < blk.aDof; ++dd) u_corr_atomic[(size_t)(blk.aOff + dd)].fetch_add(incA[dd], std::memory_order_relaxed);
        }
        if (blk.bDof > 0) {
            const VectorXr incB = MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(blk.Jb.transpose() * d);
            for (int dd = 0; dd < blk.bDof; ++dd) u_corr_atomic[(size_t)(blk.bOff + dd)].fetch_add(incB[dd], std::memory_order_relaxed);
        }
    }
}

}  // namespace

VectorXr CondensedSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::CondensedSetup);

    CondensedTopology topo = m_assembler.buildTopology();
    m_assembler.updateCompliance(topo, dt, theta);
    const VectorXr u_free = m_assembler.ufree(dt, theta);

    const int nSprings = topo.springRows;
    const int nDampers = topo.damperRows;
    const int nContacts = topo.frictionlessRows + topo.frictionalRows;
    const int nV = (int)u_free.size();

    if (m_cfg.debug_pj) {
        std::cout << "Condensed solve: nSprings = " << nSprings << ", nDampers = " << nDampers << ", nContacts = " << nContacts << ", numLambda = " << topo.numLambda << std::endl;
    }

    if (topo.numLambda == 0) return u_free;

    const VectorXr rhs = m_assembler.rhs(topo, dt, theta, u_free);
    const auto& MinvDiag = m_dyn.MinvDiag();

    // Warmstart -- identical boundary behavior to PGS (same internal sign convention, so no
    // translation needed at either boundary).
    VectorXr lambda = VectorXr::Zero(topo.numLambda);
    {
        VectorXr Lambda_g = m_dyn.Lambda_g();
        if ((int)Lambda_g.size() == nSprings) lambda.head(nSprings) = Lambda_g;
        VectorXr Lambda_gamma = m_dyn.Lambda_gamma();
        if ((int)Lambda_gamma.size() == nDampers) lambda.segment(nSprings, nDampers) = Lambda_gamma;
    }
    if (m_cfg.pj_warmstart && nContacts > 0) {
        VectorXr l_contact = VectorXr::Zero(nContacts);
        WarmstartProvider::applyWarmstart(l_contact, m_dyn, /*invertNormalSign=*/true);
        lambda.segment(nSprings + nDampers, nContacts) = l_contact;
    }

    VectorXr u_corr = VectorXr::Zero(nV);
    for (const auto& blk : topo.blocks) scatterDelta(blk, MinvDiag, lambda.segment(blk.offset, blk.dim), u_corr);

#ifdef CARDILLO_HAVE_OPENMP
    if (m_cfg.condensed_num_threads > 0) omp_set_num_threads(m_cfg.condensed_num_threads);
#endif

    sc_setup.~Scope();
    auto sc_solve = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::Condensed);

    const real_t alpha = m_cfg.pj_alpha;
    const real_t relaxation = m_cfg.pj_relaxation;

    VectorXr res_v0(nV), res_v1(nV), res_diff(nV), res_q(nV);
    auto residualNorm = [&](const VectorXr& u_prev, const VectorXr& u_cur) -> real_t {
        res_diff.noalias() = u_cur - u_prev;
        res_v0.noalias() = u_free - u_prev;
        res_v1.noalias() = u_free - u_cur;
        res_q.array() = res_v0.cwiseAbs().cwiseMax(res_v1.cwiseAbs()).array() * m_cfg.pj_tol_rel + m_cfg.pj_tol_abs;
        return (real_t)(1.0 / std::sqrt((double)nV) * res_diff.cwiseQuotient(res_q).norm());
    };

    const auto& bodyVelOffsets = m_dyn.bodyVelOffsets();
    const std::string& sweepMode = m_cfg.condensed_sweep_mode;

    cardillo::misc::Coloring coloring;
    if (sweepMode == "colored") {
        auto sc_color = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::CondensedColoring);
        const auto adjacency = buildBlockAdjacency(topo);
        coloring = cardillo::misc::colorGreedyWelshPowell((int)topo.blocks.size(), adjacency);
    }

    const bool useNewton = (m_cfg.condensed_local_solve == "newton");
    const NewtonACParams newtonParams{m_cfg.condensed_newton_max_iters, m_cfg.condensed_newton_tol};

    // Chaotic mode keeps its own persistent atomic representation of u_corr across sweeps (that's
    // the whole point -- concurrent, unsynchronized-order updates to a shared running state); it's
    // copied back into the plain `u_corr` after every sweep purely for the (single-threaded)
    // residual check and the final return value.
    std::vector<std::atomic<real_t>> u_corr_atomic;
    std::vector<int> chaoticOrder;
    std::mt19937 chaoticRng(m_cfg.condensed_chaotic_seed);
    if (sweepMode == "chaotic") {
        u_corr_atomic = std::vector<std::atomic<real_t>>(nV);
        for (int i = 0; i < nV; ++i) u_corr_atomic[(size_t)i].store(u_corr[i], std::memory_order_relaxed);
        chaoticOrder.resize(topo.blocks.size());
        std::iota(chaoticOrder.begin(), chaoticOrder.end(), 0);
        if (m_cfg.debug_pj) std::cout << "Condensed chaotic: atomic<real_t>::is_always_lock_free=" << std::atomic<real_t>::is_always_lock_free << std::endl;
    }

    VectorXr u_prev = u_corr;
    for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
        auto sc_sweep = m_dyn.timings()->scope(cardillo::misc::TimingManager::TimerId::CondensedSweep);
        if (sweepMode == "jacobi") {
            jacobiSweep(topo, bodyVelOffsets, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lambda, u_corr);
        } else if (sweepMode == "colored") {
            coloredSweep(topo, coloring, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lambda, u_corr);
        } else if (sweepMode == "chaotic") {
            chaoticSweep(topo, chaoticOrder, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, m_cfg.condensed_chaotic_reshuffle_interval, iter, chaoticRng, lambda, u_corr_atomic);
            for (int i = 0; i < nV; ++i) u_corr[i] = u_corr_atomic[(size_t)i].load(std::memory_order_relaxed);
        } else {
            gaussSeidelSweep(topo, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lambda, u_corr);
        }
        sc_sweep.~Scope();

        m_last_iters = iter + 1;
        const real_t res_norm = residualNorm(u_prev, u_corr);
        if (std::isnan(res_norm) || std::isinf(res_norm)) {
            std::cerr << "[Condensed] Divergence detected after " << iter << " iterations. Residual norm: " << res_norm << std::endl;
            throw std::runtime_error("Condensed solver diverged");
        }
        if (res_norm < (real_t)1) break;
        u_prev = u_corr;
    }

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) WarmstartProvider::storeImpulse(lambda.segment(nSprings + nDampers, nContacts), m_dyn, /*invertNormalSign=*/true);

    return u_free - u_corr;
}

}  // namespace cardillo::solver
