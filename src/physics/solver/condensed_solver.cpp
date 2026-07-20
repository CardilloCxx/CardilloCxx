#include "condensed_solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_map>

#ifdef CARDILLO_HAVE_OPENMP
#include <omp.h>
#endif

#include "../../misc/graph_coloring.hpp"
#include "local_contact_newton.hpp"
#include "warmstart.hpp"

namespace cardillo::solver {

using physics::assembly::CondensedTopology;
using physics::assembly::RowBlock;

namespace {

// Every RowBlock has dim <= 6 (1: frictionless contact, 3: frictional contact, up to 6: spring/
// damper rows) and every body has dof <= 6 (3 for a point mass, 6 for a full rigid body). Buf6 is
// a fixed-size, stack-allocated scratch buffer sized to the worst case, used with `.head(n)`
// views for whichever n <= 6 actually applies. Unlike VectorXr, constructing one never touches the
// heap -- this is what actually matters here: the per-block hot path (blockResidual/
// computeBlockUpdate/scatterDelta and their call sites below) used to allocate several VectorXr
// temporaries per block per sweep, which dominated wall-clock on scenes with many contacts (a
// malloc/free per r/lamOld/lamNew/delta/dlam, times tens of thousands of blocks, times thousands
// of sweeps). Switching those temporaries to Buf6 removes that allocation churn entirely; `Ja`/
// `Jb`/`Gii`/`GiiInv` stay dynamically-sized MatrixXXr (built once per solve() call in
// CondensedAssembler, not per sweep, so their allocation cost is amortized and not the bottleneck).
using Buf6 = Vector6r;

// Same clip as ProjectedGaussSeidel's project()/project_all(): normal impulse (index 0) clamped to
// <= 0, tangential (indices 1,2) scaled toward the Coulomb disk of radius -mu*lambda[0]. No-op for
// Spring/Damper blocks (equality constraints, unconstrained). Templated so it works identically on
// a VectorXr::segment() view or a Buf6::head() view -- both are Eigen::MatrixBase-derived.
template <typename Derived>
void projectBlock(RowBlock::Kind kind, real_t mu, Eigen::MatrixBase<Derived>& lam) {
    if (kind == RowBlock::Kind::Spring || kind == RowBlock::Kind::Damper) return;

    lam[0] = std::min(lam[0], (real_t)0);
    if (kind == RowBlock::Kind::ContactFrictionless) return;

    const real_t tnorm = std::sqrt(lam[1] * lam[1] + lam[2] * lam[2]);
    if (tnorm <= (real_t)0) return;
    const real_t s = std::min((real_t)1, (-mu * lam[0]) / tnorm);
    lam[1] *= s;
    lam[2] *= s;
}

// Writes this block's local residual (rhs_block - Ja*u_corr_A - Jb*u_corr_B -
// complianceDiag*lambda_block) into out.head(blk.dim). No heap allocation: `out` is caller-owned
// stack storage, written into via `.noalias()` rather than returned by value.
void blockResidual(const RowBlock& blk, const VectorXr& rhs, const VectorXr& u_corr, const VectorXr& lambda, Buf6& out) {
    auto o = out.head(blk.dim);
    o = rhs.segment(blk.offset, blk.dim);
    if (blk.aDof > 0) o.noalias() -= blk.Ja * u_corr.segment(blk.aOff, blk.aDof);
    if (blk.bDof > 0) o.noalias() -= blk.Jb * u_corr.segment(blk.bOff, blk.bDof);
    o.noalias() -= blk.complianceDiag.cwiseProduct(lambda.segment(blk.offset, blk.dim));
}

// nullptr for every body without an active implicit-gyroscopic override -- the overwhelming common
// case, so every call site below keeps its original plain-diagonal-MinvDiag expression unchanged in
// that branch. See CondensedTopology::gyroMinvBlocks.
const Matrix66r* gyroBlockFor(int bodyIndex, const std::unordered_map<int, Matrix66r>& gyroBlocks) {
    if (bodyIndex < 0 || gyroBlocks.empty()) return nullptr;
    auto it = gyroBlocks.find(bodyIndex);
    return it == gyroBlocks.end() ? nullptr : &it->second;
}

// Applies one block's impulse delta (dlambda.head(blk.dim)) to u_corr in place: u_corr += Minv *
// J^T * dlambda. `tmp` is caller-provided scratch for the intermediate J^T*dlambda product (also
// Buf6, also heap-free). `gyroBlocks` is topo.gyroMinvBlocks -- empty for every scene without an
// active implicit-gyroscopic body, in which case every branch below reduces to exactly the original
// expression (see gyroBlockFor()).
void scatterDelta(const RowBlock& blk, const VectorXr& MinvDiag, const std::unordered_map<int, Matrix66r>& gyroBlocks, const Buf6& dlambda, VectorXr& u_corr, Buf6& tmp) {
    if (blk.aDof > 0) {
        tmp.head(blk.aDof).noalias() = blk.Ja.transpose() * dlambda.head(blk.dim);
        if (const auto* gA = gyroBlockFor(blk.bodyIndexA, gyroBlocks))
            u_corr.segment(blk.aOff, blk.aDof).noalias() += gA->topLeftCorner(blk.aDof, blk.aDof) * tmp.head(blk.aDof);
        else
            u_corr.segment(blk.aOff, blk.aDof).noalias() += MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(tmp.head(blk.aDof));
    }
    if (blk.bDof > 0) {
        tmp.head(blk.bDof).noalias() = blk.Jb.transpose() * dlambda.head(blk.dim);
        if (const auto* gB = gyroBlockFor(blk.bodyIndexB, gyroBlocks))
            u_corr.segment(blk.bOff, blk.bDof).noalias() += gB->topLeftCorner(blk.bDof, blk.bDof) * tmp.head(blk.bDof);
        else
            u_corr.segment(blk.bOff, blk.bDof).noalias() += MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(tmp.head(blk.bDof));
    }
}

// condensed.true_schur=true: exactly zeroes every bilateral (spring+damper) row's residual in one
// shot, given the CURRENT contact contribution to u_corr, via the precomputed block-sparse LDLT of
// Sbb (see CondensedAssembler::buildBilateralFactorization()). This replaces those rows'
// participation in the Gauss-Seidel/Jacobi/colored sweep entirely -- bilateral rows are pure
// equality constraints (projectBlock() already no-ops for them), so unlike the projected contact
// rows there is no iteration needed once their joint linear system is solved exactly. A no-op when
// there are no bilateral blocks (e.g. domino), by construction.
void exactBilateralStep(const CondensedTopology& topo, const misc::BlockSparseLDLT& ldlt, const VectorXr& rhs, const VectorXr& MinvDiag, VectorXr& lambda, VectorXr& u_corr) {
    const int n = topo.numBilateralBlocks;
    if (n == 0) return;

    VectorXr r(ldlt.totalDim());
    Buf6 rb;
    for (int i = 0; i < n; ++i) {
        const auto& blk = topo.blocks[(size_t)i];
        blockResidual(blk, rhs, u_corr, lambda, rb);
        r.segment(blk.offset, blk.dim) = rb.head(blk.dim);
    }

    const VectorXr delta = ldlt.solve(r);

    Buf6 d, tmp;
    for (int i = 0; i < n; ++i) {
        const auto& blk = topo.blocks[(size_t)i];
        d.head(blk.dim) = delta.segment(blk.offset, blk.dim);
        scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, d, u_corr, tmp);
        lambda.segment(blk.offset, blk.dim) += d.head(blk.dim);
    }
}

bool isContactKind(RowBlock::Kind kind) { return kind == RowBlock::Kind::ContactFrictionless || kind == RowBlock::Kind::ContactFrictional; }

// Shared local-solve step used by every sweep driver: for frictional contact blocks when
// `useNewton` is set, tries the guarded Alart-Curnier semismooth Newton solve first (see
// local_contact_newton.hpp) -- it either converges to a locally-exact solution in a handful of
// iterations, or reports failure (singular Jacobian, non-finite step, non-decreasing residual, or
// no convergence within its iteration cap), in which case we fall through to the always-available
// linear step + projection below. Writes the new (projected/converged) impulse into
// out.head(blk.dim); caller computes the delta.
// `linearOnly` (default false, so every existing call site is unaffected): skips both the Newton
// branch and the final projectBlock() call, making this update a pure linear relaxation step
// (o = relaxation*alpha*GiiInv*r + lamOld, nothing else). Used only by estimateSpectralRadius()
// below, which needs the *linear* part of the sweep operator in isolation -- Newton/projection are
// exactly the nonlinearities that make a rigorous spectral-radius argument inapplicable in the
// first place (see that function's own comment).
void computeBlockUpdate(const RowBlock& blk, const Buf6& r, real_t relaxation, real_t alpha, const Buf6& lamOld, bool useNewton, const NewtonACParams& newtonParams, Buf6& out,
                         bool linearOnly = false) {
    if (useNewton && !linearOnly && blk.kind == RowBlock::Kind::ContactFrictional) {
        Vector3r lam3 = lamOld.head<3>();
        const Vector3r r3 = r.head<3>();
        if (solveContactBlockNewtonAC(blk.Gii, r3, blk.mu, lam3, newtonParams)) {
            out.head<3>() = lam3;
            return;
        }
    }

    auto o = out.head(blk.dim);
    o.noalias() = relaxation * (blk.GiiInv * r.head(blk.dim));
    if (isContactKind(blk.kind)) o *= alpha;
    o += lamOld.head(blk.dim);
    if (!linearOnly) projectBlock(blk.kind, blk.mu, o);
}

// Sequential block-Gauss-Seidel sweep: propagates u_corr immediately after each block, matching
// ProjectedGaussSeidelSolver's own loop structure exactly, just expressed over RowBlocks.
// `startBlock` skips the bilateral (spring+damper) prefix when condensed.true_schur handles those
// exactly instead (0 -- i.e. every block -- otherwise); contact blocks are always a contiguous
// suffix of topo.blocks, so a single start index is enough, no separate index list needed.
void gaussSeidelSweep(const CondensedTopology& topo, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton, const NewtonACParams& newtonParams,
                       VectorXr& lambda, VectorXr& u_corr, int startBlock = 0, bool linearOnly = false) {
    Buf6 r, lamOld, lamNew, delta, tmp;
    for (size_t bi = (size_t)startBlock; bi < topo.blocks.size(); ++bi) {
        const auto& blk = topo.blocks[bi];
        blockResidual(blk, rhs, u_corr, lambda, r);
        lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
        computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew, linearOnly);

        delta.head(blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
        scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, delta, u_corr, tmp);
        lambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim);
    }
}

// Builds the block adjacency graph for coloring: two blocks are adjacent iff they share a body
// (i.e. both appear in that body's incident-block list from CondensedTopology::blocksOfBody).
// Deliberately NOT the dense Schur-complement/Delassus sparsity pattern (which for a compliant
// chain would be fully dense and collapse coloring into one color) -- this is the same
// shared-body adjacency Vivace/Bullet-style contact coloring uses.
//
// `startBlock` excludes the bilateral (spring+damper) prefix when condensed.true_schur is active
// (those rows are handled by an exact solve, not colored/swept at all): the returned adjacency is
// LOCAL (0-based over just blocks[startBlock:]), matching what colorGreedyWelshPowell() expects --
// callers must add startBlock back to every index in the resulting color classes before using them
// against topo.blocks.
std::vector<std::vector<int>> buildBlockAdjacency(const CondensedTopology& topo, int startBlock = 0) {
    std::vector<std::vector<int>> adj(topo.blocks.size() - (size_t)startBlock);
    for (const auto& incident : topo.blocksOfBody) {
        std::vector<int> local;
        local.reserve(incident.size());
        for (int idx : incident) {
            if (idx >= startBlock) local.push_back(idx - startBlock);
        }
        for (size_t i = 0; i < local.size(); ++i) {
            for (size_t j = i + 1; j < local.size(); ++j) {
                adj[(size_t)local[i]].push_back(local[j]);
                adj[(size_t)local[j]].push_back(local[i]);
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
void coloredSweep(const CondensedTopology& topo, const misc::Coloring& coloring, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                   const NewtonACParams& newtonParams, VectorXr& lambda, VectorXr& u_corr, bool linearOnly = false) {
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
            // Declared inside the loop body: each thread gets its own stack-local (and heap-free)
            // instance automatically, no `private()` clause needed.
            Buf6 r, lamOld, lamNew, delta, tmp;
            const int idx = cls[(size_t)k];
            const auto& blk = topo.blocks[(size_t)idx];
            blockResidual(blk, rhs, u_corr, lambda, r);
            lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
            computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew, linearOnly);

            // Safe: no two blocks in `cls` share a body (that's the coloring invariant), so these
            // writes land in disjoint slices of u_corr/lambda across threads -- no race, no atomics.
            delta.head(blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, delta, u_corr, tmp);
            lambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim);
        }
        // implicit barrier at the end of `omp for` here, before the next color starts
    }
}

// Jacobi sweep: two passes, no shared-state races even when parallelized (see design note in
// condensed_solver.hpp). Pass 1 updates each block's own impulse against a frozen snapshot of
// u_corr from the previous sweep. Pass 2 has each BODY gather corrections from only its own
// incident blocks -- never a block scattering into shared state -- so the reduction needs no
// synchronization at all, not even atomics.
// `startBlock` skips the bilateral prefix under condensed.true_schur, same convention as
// gaussSeidelSweep/buildBlockAdjacency above. Pass 2's gather is left iterating every body's full
// incident list unchanged: dlambda for any skipped (bilateral) block stays at its zero-init value
// from `VectorXr::Zero(lambda.size())` below, since pass 1 never touches that segment, so gathering
// it contributes nothing -- restricting pass 2's range too would save a little work but isn't
// needed for correctness.
void jacobiSweep(const CondensedTopology& topo, const std::vector<int>& bodyVelOffsets, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                  const NewtonACParams& newtonParams, VectorXr& lambda, VectorXr& u_corr, int startBlock = 0, bool linearOnly = false) {
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
        for (int i = startBlock; i < numBlocks; ++i) {
            Buf6 r, lamOld, lamNew;
            const auto& blk = topo.blocks[i];
            blockResidual(blk, rhs, u_read, lambda, r);
            lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
            computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew, linearOnly);

            dlambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
            lambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim);
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

            Buf6 acc = Buf6::Zero(), d;
            for (int blockIdx : incident) {
                const auto& blk = topo.blocks[blockIdx];
                d.head(blk.dim) = dlambda.segment(blk.offset, blk.dim);
                if (blk.bodyIndexA == b) acc.head(dof).noalias() += blk.Ja.transpose() * d.head(blk.dim);
                if (blk.bodyIndexB == b) acc.head(dof).noalias() += blk.Jb.transpose() * d.head(blk.dim);
            }
            if (const auto* gB = gyroBlockFor(b, topo.gyroMinvBlocks))
                u_write.segment(off, dof) = u_read.segment(off, dof) + gB->topLeftCorner(dof, dof) * acc.head(dof);
            else
                u_write.segment(off, dof) = u_read.segment(off, dof) + MinvDiag.segment(off, dof).cwiseProduct(acc.head(dof));
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
        Buf6 uA = Buf6::Zero(), uB = Buf6::Zero(), r, lamOld, lamNew, d, inc;
        const int idx = order[(size_t)k];
        const auto& blk = topo.blocks[(size_t)idx];

        for (int dd = 0; dd < blk.aDof; ++dd) uA[dd] = u_corr_atomic[(size_t)(blk.aOff + dd)].load(std::memory_order_relaxed);
        for (int dd = 0; dd < blk.bDof; ++dd) uB[dd] = u_corr_atomic[(size_t)(blk.bOff + dd)].load(std::memory_order_relaxed);

        auto ro = r.head(blk.dim);
        ro = rhs.segment(blk.offset, blk.dim);
        if (blk.aDof > 0) ro.noalias() -= blk.Ja * uA.head(blk.aDof);
        if (blk.bDof > 0) ro.noalias() -= blk.Jb * uB.head(blk.bDof);
        ro.noalias() -= blk.complianceDiag.cwiseProduct(lambda.segment(blk.offset, blk.dim));

        lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
        computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew);
        lambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim);  // safe: `order` is a permutation this sweep

        d.head(blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
        if (blk.aDof > 0) {
            inc.head(blk.aDof).noalias() = blk.Ja.transpose() * d.head(blk.dim);
            if (const auto* gA = gyroBlockFor(blk.bodyIndexA, topo.gyroMinvBlocks))
                inc.head(blk.aDof) = gA->topLeftCorner(blk.aDof, blk.aDof) * inc.head(blk.aDof);
            else
                inc.head(blk.aDof) = MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(inc.head(blk.aDof));
            for (int dd = 0; dd < blk.aDof; ++dd) u_corr_atomic[(size_t)(blk.aOff + dd)].fetch_add(inc[dd], std::memory_order_relaxed);
        }
        if (blk.bDof > 0) {
            inc.head(blk.bDof).noalias() = blk.Jb.transpose() * d.head(blk.dim);
            if (const auto* gB = gyroBlockFor(blk.bodyIndexB, topo.gyroMinvBlocks))
                inc.head(blk.bDof) = gB->topLeftCorner(blk.bDof, blk.bDof) * inc.head(blk.bDof);
            else
                inc.head(blk.bDof) = MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(inc.head(blk.bDof));
            for (int dd = 0; dd < blk.bDof; ++dd) u_corr_atomic[(size_t)(blk.bOff + dd)].fetch_add(inc[dd], std::memory_order_relaxed);
        }
    }
}

// Estimates the spectral radius of the LINEARIZED sweep operator (i.e. with projectBlock()/Newton
// disabled -- exactly the nonlinearities that make a rigorous spectral-radius argument
// inapplicable in the first place, same caveat as ProjectedJacobiSolver's own estimator), for
// Chebyshev semi-iterative acceleration (condensed.pj_chebyshev... see pj.chebyshev). Mirrors
// ProjectedJacobiSolver's own estimateSpectralRadius(), adapted to condensed's matrix-free
// block-sweep representation instead of a global sparse system: rather than needing an explicit
// linear operator/solveS(), this runs the SAME sweep function the real iteration uses
// (gaussSeidelSweep/jacobiSweep/coloredSweep, plus exactBilateralStep under true_schur) with
// projection/Newton disabled (linearOnly=true) and a ZERO rhs. With no projection and no rhs, one
// sweep step is an EXACT linear map G(d) = A*d (no affine offset -- the rhs-dependent constant
// term vanishes), so plain power iteration on a probe vector converges to |A|'s dominant
// eigenvalue magnitude, exactly analogous to PJ's own power iteration on its (implicit,
// matrix-based) linearized operator. Uses a deterministic all-ones probe vector (matching PJ's own
// choice) rather than a randomized one, so this doesn't introduce a new nondeterminism source into
// the solver.
real_t estimateSpectralRadius(const CondensedTopology& topo, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, const std::string& sweepMode, const misc::Coloring& coloring,
                               const std::vector<int>& bodyVelOffsets, int startBlock, bool useTrueSchur, const misc::BlockSparseLDLT& bilateralLdlt, int nV) {
    constexpr int kPowerIterations = 8;
    const int numLambda = topo.numLambda;
    if (numLambda == 0) return (real_t)0;

    const VectorXr zeroRhs = VectorXr::Zero(numLambda);
    const NewtonACParams unusedNewtonParams{};  // never consulted: linearOnly=true skips the Newton branch entirely

    VectorXr d_lambda = VectorXr::Ones(numLambda).normalized();
    VectorXr d_u_corr = VectorXr::Zero(nV);
    {
        Buf6 lam0, tmp0;
        for (const auto& blk : topo.blocks) {
            lam0.head(blk.dim) = d_lambda.segment(blk.offset, blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, lam0, d_u_corr, tmp0);
        }
    }

    real_t rho = 0;
    for (int it = 0; it < kPowerIterations; ++it) {
        if (useTrueSchur) exactBilateralStep(topo, bilateralLdlt, zeroRhs, MinvDiag, d_lambda, d_u_corr);
        if (sweepMode == "jacobi") {
            jacobiSweep(topo, bodyVelOffsets, zeroRhs, MinvDiag, relaxation, alpha, /*useNewton=*/false, unusedNewtonParams, d_lambda, d_u_corr, startBlock, /*linearOnly=*/true);
        } else if (sweepMode == "colored") {
            coloredSweep(topo, coloring, zeroRhs, MinvDiag, relaxation, alpha, false, unusedNewtonParams, d_lambda, d_u_corr, /*linearOnly=*/true);
        } else {
            gaussSeidelSweep(topo, zeroRhs, MinvDiag, relaxation, alpha, false, unusedNewtonParams, d_lambda, d_u_corr, startBlock, /*linearOnly=*/true);
        }

        const real_t n = d_lambda.norm();
        // A genuine bug found by tracing this per-iteration on slinky's first (contact-free) step: when
        // condensed.true_schur exactly solves the ENTIRE probe (e.g. no contacts at all, so the
        // whole linear-only sweep is just exactBilateralStep on a homogeneous system), `n` collapses
        // to floating-point noise (~1e-13, confirmed by direct trace), not a real small eigenvalue.
        // The old `1e-300` threshold is nowhere near the actual FP noise floor for this
        // O(1)-scale computation, so it let that noise get renormalized back to unit length and
        // amplified ~10^12x every remaining iteration -- converging to a spurious rho close to 1
        // (an artifact of amplified rounding error, not a real eigenvalue) instead of correctly
        // reporting "this probe found nothing to accelerate". 1e-8 sits safely above realistic FP
        // noise for this scale and well below any genuine eigenvalue magnitude worth resolving.
        if (!(n > (real_t)1e-8)) return (real_t)0;
        rho = n;
        d_lambda /= n;
        d_u_corr /= n;
    }
    // Deliberately NOT clamped here -- callers need the raw value for different purposes (the
    // Chebyshev omega recurrence needs it clamped strictly below 1; a caller deriving a safe `alpha`
    // from this estimate needs the true, possibly-much-larger-than-1 magnitude to know how far past
    // stability the current alpha/relaxation choice actually is). Clamp at each call site instead.
    if (getenv("CARDILLO_DEBUG_RHO_RAW")) std::cout << "[Condensed] raw (unclamped) spectral radius estimate: " << rho << std::endl;
    return rho;
}

}  // namespace

VectorXr CondensedSolver::solve(real_t dt, real_t theta) {
    auto sc_setup = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);

    CondensedTopology topo = m_assembler.buildTopology(dt);
    m_assembler.updateCompliance(topo, dt, theta);
    const VectorXr u_free = m_assembler.ufree(dt, theta, topo);

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
    {
        Buf6 lam0, tmp0;
        for (const auto& blk : topo.blocks) {
            lam0.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, lam0, u_corr, tmp0);
        }
    }

#ifdef CARDILLO_HAVE_OPENMP
    if (m_cfg.condensed_num_threads > 0) omp_set_num_threads(m_cfg.condensed_num_threads);
#endif

    sc_setup.~Scope();
    auto sc_solve = m_dyn.timings()->scope(misc::TimingManager::TimerId::Condensed);

    real_t alpha = m_cfg.pj_alpha;
    const real_t relaxation = m_cfg.pj_relaxation;

    VectorXr res_v0(nV), res_v1(nV), res_diff(nV), res_q(nV);
    auto residualNorm = [&](const VectorXr& u_prev, const VectorXr& u_cur) -> real_t {
        res_diff.noalias() = u_cur - u_prev;
        res_v0.noalias() = u_free - u_prev;
        res_v1.noalias() = u_free - u_cur;
        res_q.array() = res_v0.cwiseAbs().cwiseMax(res_v1.cwiseAbs()).array() * m_cfg.pj_tol_rel + m_cfg.pj_tol_abs;
        return (real_t)(1.0 / std::sqrt((double)nV) * res_diff.cwiseQuotient(res_q).norm());
    };
    auto checkResidual = [&](const VectorXr& uA, const VectorXr& uB, int iter) -> real_t {
        const real_t res_norm = residualNorm(uA, uB);
        if (std::isnan(res_norm) || std::isinf(res_norm)) {
            std::cerr << "[Condensed] Divergence detected after " << iter << " iterations. Residual norm: " << res_norm << std::endl;
            throw std::runtime_error("Condensed solver diverged");
        }
        return res_norm;
    };

    const auto& bodyVelOffsets = m_dyn.bodyVelOffsets();
    const std::string& sweepMode = m_cfg.condensed_sweep_mode;

    // condensed.true_schur only applies to gauss_seidel/jacobi/colored -- chaotic keeps its own
    // atomic u_corr representation and full block range unchanged (see the doSweep comment below
    // for why chaotic is excluded from the generic wrapper in the first place; combining chaotic's
    // deliberately-stale atomic updates with an exact joint bilateral solve raises questions not
    // tackled here). `startBlock` is where contact blocks begin -- 0 (i.e. every block still swept
    // together) unless true_schur is on, in which case it's topo.numBilateralBlocks.
    const bool useTrueSchur = m_cfg.condensed_true_schur && sweepMode != "chaotic";
    const int startBlock = useTrueSchur ? topo.numBilateralBlocks : 0;
    if (m_cfg.condensed_true_schur && sweepMode == "chaotic" && m_cfg.debug_pj) {
        std::cout << "[Condensed] condensed.true_schur has no effect on sweep_mode=chaotic; ignoring." << std::endl;
    }

    misc::BlockSparseLDLT bilateralLdlt;
    if (useTrueSchur && topo.numBilateralBlocks > 0) {
        auto sc_schur = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);
        bilateralLdlt = m_assembler.buildBilateralFactorization(topo);
    }

    misc::Coloring coloring;
    if (sweepMode == "colored") {
        auto sc_color = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedColoring);
        const auto adjacency = buildBlockAdjacency(topo, startBlock);
        coloring = misc::colorGreedyWelshPowell((int)topo.blocks.size() - startBlock, adjacency);
        for (auto& cls : coloring.colorClasses)
            for (auto& idx : cls) idx += startBlock;
    }

    const bool useNewton = (m_cfg.condensed_local_solve == "newton");
    const NewtonACParams newtonParams{m_cfg.condensed_newton_max_iters, m_cfg.condensed_newton_tol,
                                       (m_cfg.condensed_newton_rho_strategy == "full") ? NewtonRhoStrategy::FullSpectral : NewtonRhoStrategy::Split};

    // condensed.auto_alpha (experimental, default off -- see config.hpp): derive `alpha` for THIS
    // step from the same power-iteration estimator pj.chebyshev uses, instead of the fixed
    // hand-tuned pj.alpha. Measures the raw spectral radius at alpha=1 (rho1), then picks alpha so
    // the resulting operator's radius lands near condensed_auto_alpha_target_rho, using the
    // empirically-confirmed linear relationship rho(alpha) ~= alpha*(rho1+1) - 1 for alpha large
    // enough to be dominated by the largest eigenvalue (see CONDENSED_SOLVER_REPORT.md). Overrides
    // `alpha` for every use below, including the Chebyshev estimate itself if also enabled.
    if (m_cfg.condensed_auto_alpha) {
        auto sc_autoalpha = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);
        const real_t rho1 = estimateSpectralRadius(topo, MinvDiag, relaxation, (real_t)1.0, sweepMode, coloring, bodyVelOffsets, startBlock, useTrueSchur, bilateralLdlt, nV);
        if (rho1 > 0) {
            const real_t targetRho = m_cfg.condensed_auto_alpha_target_rho;
            alpha = std::clamp((targetRho + (real_t)1.0) / (rho1 + (real_t)1.0), (real_t)1e-4, (real_t)1.0);
        }
        if (m_cfg.debug_pj) std::cout << "[Condensed] auto_alpha: raw rho(alpha=1)=" << rho1 << ", derived alpha=" << alpha << std::endl;
    }

    // condensed.pj_chebyshev reuses PJ's own config key (pj.chebyshev) and precedence rule --
    // pj_nesterov takes priority if both are set (matching ProjectedJacobiSolver::solve()'s
    // `useChebyshev = pj_chebyshev && !pj_nesterov`) -- rather than inventing a condensed-specific
    // key or a new precedence order. Excluded on sweep_mode=chaotic, same reasoning as true_schur
    // above: chaotic's atomic u_corr representation isn't compatible with extrapolating a plain
    // VectorXr state between calls.
    const bool useChebyshev = m_cfg.pj_chebyshev && !m_cfg.pj_nesterov && sweepMode != "chaotic";
    if (m_cfg.pj_chebyshev && !m_cfg.pj_nesterov && sweepMode == "chaotic" && m_cfg.debug_pj) {
        std::cout << "[Condensed] pj.chebyshev has no effect on sweep_mode=chaotic; ignoring." << std::endl;
    }
    real_t chebyshevRho = 0;
    if (useChebyshev && !m_cfg.condensed_chebyshev_adaptive_rho) {
        // Skipped entirely when condensed_chebyshev_adaptive_rho is set: that path derives rho from
        // a warmup of real sweeps instead (see the Chebyshev branch below), making this upfront
        // linearized power-iteration estimate redundant -- no point paying its cost just to
        // immediately overwrite the result.
        //
        // Chebyshev's spectral-radius estimate is only rigorously meaningful when the linearized
        // sweep operator has a REAL eigenvalue spectrum -- guaranteed (via similarity to a
        // symmetric matrix) when Minv is symmetric, i.e. no active implicit-gyroscopic override
        // anywhere in the system (see docs/chapters/solvers/condensed.rst's Chebyshev section).
        // A non-symmetric Minv can give the operator complex eigenvalues, which the classical
        // 3-term Chebyshev recurrence isn't designed for -- not a guaranteed failure (the same
        // "heuristic on the full nonlinear problem" caveat already applies regardless), but worth
        // surfacing rather than silently trusting the estimate. Warn only (debug_pj-gated, like
        // every other diagnostic here) rather than auto-disabling: unlike true_schur+chaotic
        // (a hard structural incompatibility), this is a theoretical-grounding caveat with no
        // empirical evidence either way yet, so changing behavior here isn't justified.
        if (m_cfg.debug_pj && !topo.gyroMinvBlocks.empty()) {
            std::cout << "[Condensed] Warning: pj.chebyshev's spectral-radius estimate assumes a real "
                         "eigenvalue spectrum (via symmetric Minv); an implicit-gyroscopic body is active, "
                         "making Minv non-symmetric -- the estimate may not be reliable here."
                      << std::endl;
        }
        auto sc_cheb = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);
        const real_t rawRho = estimateSpectralRadius(topo, MinvDiag, relaxation, alpha, sweepMode, coloring, bodyVelOffsets, startBlock, useTrueSchur, bilateralLdlt, nV);
        // Clamp strictly below 1: the omega recurrence divides by (1 - rho^2/4*omega), only
        // well-defined/stable for rho < 1 (a converging basic iteration to begin with) -- same
        // clamp PJ's own estimateSpectralRadius() uses. estimateSpectralRadius() itself returns the
        // raw, unclamped value now (see its own comment) -- clamping is this call site's job.
        chebyshevRho = std::min(rawRho, (real_t)0.995);
        if (m_cfg.debug_pj) std::cout << "[Condensed] chebyshev spectral radius estimate: " << chebyshevRho << std::endl;
    }

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

    // `doSweep` is generic over sweep mode so Nesterov acceleration (below) can wrap any of
    // gauss_seidel/jacobi/colored without knowing their internals -- it only ever touches the
    // outer (lambda, u_corr) state between calls, exactly like ProjectedGaussSeidelSolver's own
    // Nesterov branch does around pgs_sweep(). `chaotic` is deliberately excluded: its state lives
    // in `u_corr_atomic`, not a plain VectorXr, and momentum-extrapolating a value that's already
    // being intentionally kept stale by design is a separate question not tackled here.
    //
    // Under condensed.true_schur, every call first exactly solves the bilateral (spring+damper)
    // rows given the sweep's current contact contribution to u_corr, then runs the usual sweep
    // restricted to the contact blocks (startBlock) -- the bilateral rows never enter the sweep
    // itself. `useTrueSchur` is false (and startBlock is 0) whenever there's nothing to gain from
    // this (no bilateral rows, or true_schur not requested), making doSweep identical to before.
    auto doSweep = [&](VectorXr& lam, VectorXr& ucorr) {
        if (useTrueSchur) exactBilateralStep(topo, bilateralLdlt, rhs, MinvDiag, lam, ucorr);
        if (sweepMode == "jacobi") {
            jacobiSweep(topo, bodyVelOffsets, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lam, ucorr, startBlock);
        } else if (sweepMode == "colored") {
            coloredSweep(topo, coloring, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lam, ucorr);
        } else {
            gaussSeidelSweep(topo, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lam, ucorr, startBlock);
        }
    };

    if (sweepMode == "chaotic") {
        VectorXr u_prev = u_corr;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            chaoticSweep(topo, chaoticOrder, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, m_cfg.condensed_chaotic_reshuffle_interval, iter, chaoticRng, lambda, u_corr_atomic);
            for (int i = 0; i < nV; ++i) u_corr[i] = u_corr_atomic[(size_t)i].load(std::memory_order_relaxed);
            sc_sweep.~Scope();

            m_last_iters = iter + 1;
            if (checkResidual(u_prev, u_corr, iter) < (real_t)1) break;
            u_prev = u_corr;
        }
    } else if (!m_cfg.pj_nesterov && !useChebyshev) {
        VectorXr u_prev = u_corr;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            doSweep(lambda, u_corr);
            sc_sweep.~Scope();

            m_last_iters = iter + 1;
            if (checkResidual(u_prev, u_corr, iter) < (real_t)1) break;
            u_prev = u_corr;
        }
    } else if (m_cfg.pj_nesterov) {
        // Nesterov (FISTA-style) acceleration, ported verbatim from
        // ProjectedGaussSeidelSolver::solve()'s Nesterov branch: momentum-extrapolate between
        // sweeps, with adaptive restart when the residual grows, the momentum coefficient gets too
        // large, or the update direction reverses. Reuses the same shared config keys PJ/PGS
        // already use (pj.nesterov, pj.nesterov_beta_threshold, pj.nesterov_restart_limit) rather
        // than inventing condensed-specific ones.
        VectorXr lambda_k = lambda, u_k = u_corr, lambda_y = lambda, u_y = u_corr;
        double thk = 1.0;
        real_t err_prev = std::numeric_limits<real_t>::infinity();
        int restarts = 0;
        bool momentum_disabled = false;

        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            lambda = lambda_y;
            u_corr = u_y;
            doSweep(lambda, u_corr);
            sc_sweep.~Scope();

            const VectorXr lambda_k1 = lambda;
            const VectorXr u_k1 = u_corr;
            // Exactly one doSweep() happened this iteration -- same as the plain loop above and
            // PJ's own nesterov_loop() (which uses iter+1) -- so this must match. (Previously
            // iter+2, inherited verbatim from PGS's Nesterov branch, which had the same stray
            // off-by-one with no second sweep anywhere to justify it -- see that fix. Purely
            // cosmetic: only ever read by the progress-bar/log "iters=" display, never by any
            // convergence or control-flow decision.)
            m_last_iters = iter + 1;

            // Non-throwing here, deliberately -- matching PJ's own nesterov_loop() (see
            // projected_jacobi.cpp), which never throws on a non-finite residual either: a NaN/Inf
            // `err` is one of the restart triggers below, not an immediate hard failure. Throwing
            // here (as an earlier version of this code did via checkResidual()) made that restart
            // branch permanently dead code -- momentum could never recover from an overshoot because
            // the exception fired before the restart logic ran at all. This was the actual cause of
            // the jacobi+Nesterov divergence on slinky documented in condensed.rst/the solver report:
            // not an unrecoverable instability, just a self-inflicted inability to ever try recovering.
            const real_t err = residualNorm(u_y, u_k1);
            if (err <= (real_t)1) break;

            const double thk1 = 0.5 * (1.0 + std::sqrt(4.0 * thk * thk + 1.0));
            const double betak1 = (thk - 1.0) / thk1;

            bool restart = false;
            if (!std::isfinite((double)err)) {
                restart = true;
            } else if (err_prev < std::numeric_limits<real_t>::infinity() && err > (real_t)1.05 * err_prev) {
                restart = true;
            } else if (betak1 > (double)m_cfg.pj_nesterov_beta_threshold) {
                restart = true;
            } else if (!std::isfinite(betak1) || betak1 < 0.0 || betak1 > 1.0) {
                restart = true;
            } else if ((u_y - u_k1).dot(u_k1 - u_k) > 0.0) {
                restart = true;
            }

            // A restart only ever resets to lambda_k1/u_k1 -- the plain (unaccelerated) sweep
            // output. If that itself is already non-finite (the corruption happened one iteration
            // earlier, inside the extrapolation below, and propagated through doSweep), no restart
            // can recover it: momentum is permanently disabled and the loop's exit condition is a
            // hard failure at the bottom, exactly like the plain (non-Nesterov) path's checkResidual.
            if (!lambda_k1.allFinite() || !u_k1.allFinite()) {
                momentum_disabled = true;
                restarts = m_cfg.pj_nesterov_restart_limit;
                checkResidual(u_y, u_k1, iter);  // throws: nothing left to fall back to
            }

            if (restart) {
                lambda_y = lambda_k1;
                u_y = u_k1;
                lambda_k = lambda_k1;
                thk = 1.0;
                if (++restarts >= m_cfg.pj_nesterov_restart_limit) momentum_disabled = true;
            } else if (momentum_disabled) {
                lambda_y = lambda_k1;
                u_y = u_k1;
                thk = 1.0;
            } else {
                // Guard the extrapolation itself, not just its aftermath: betak1 in [0,1] bounds the
                // *coefficient*, not the extrapolated magnitude -- (lambda_k1 - lambda_k) can still be
                // large enough to overflow. Catching that here, before it is ever fed into the next
                // doSweep(), is what actually prevents the corruption from entering the state at all;
                // the `allFinite()` check above is only a last-resort backstop for anything that slips
                // past this guard some other way.
                VectorXr trial_lambda_y = lambda_k1 + (real_t)betak1 * (lambda_k1 - lambda_k);
                VectorXr trial_u_y = u_k1 + (real_t)betak1 * (u_k1 - u_k);
                if (trial_lambda_y.allFinite() && trial_u_y.allFinite()) {
                    lambda_y = std::move(trial_lambda_y);
                    u_y = std::move(trial_u_y);
                    thk = thk1;
                } else {
                    lambda_y = lambda_k1;
                    u_y = u_k1;
                    thk = thk1;
                }
            }
            lambda_k = lambda_k1;
            u_k = u_k1;
            err_prev = err;
        }
        // `lambda`/`u_corr` already hold the latest post-sweep values from inside the loop above
        // (doSweep mutates them by reference) -- no final reassignment needed, matching
        // ProjectedGaussSeidelSolver's own Nesterov branch, which relies on the same fact.
    } else {
        // useChebyshev is the only remaining possibility here (the first branch above excluded
        // !pj_nesterov && !useChebyshev, the second excluded pj_nesterov).
        //
        // Chebyshev semi-iterative acceleration, ported from ProjectedJacobiSolver's
        // chebyshev_loop(): unlike Nesterov, a FIXED extrapolation schedule (no residual-dependent
        // momentum/restart logic) driven entirely by chebyshevRho, the linearized sweep operator's
        // spectral radius estimated once above via estimateSpectralRadius(). Mutually exclusive
        // with Nesterov (useChebyshev already required !pj_nesterov, matching
        // ProjectedJacobiSolver::solve()'s own precedence: pj_nesterov > pj_chebyshev).
        real_t rho = chebyshevRho;
        bool warmupConverged = false;
        int warmupItersUsed = 0;
        if (m_cfg.condensed_chebyshev_adaptive_rho) {
            // Adaptive rho (the classical "adaptive SOR/Chebyshev" technique -- e.g. Manteuffel
            // 1977/Ashby 1985's power-method-based adaptive parameter estimation): rather than a
            // separate linearized, zero-rhs power-iteration pre-pass (measured earlier to be
            // overly pessimistic for contact-rich scenes -- it ignores the stabilizing effect of
            // the friction-cone projection entirely, treating every contact as permanently,
            // fully coupled), run a few REAL (unaccelerated, actually-projected) sweeps first and
            // estimate rho from the OBSERVED residual ratio: for a linearly-converging iteration,
            // err_k/err_{k-1} -> rho asymptotically. These sweeps are real solve work, not thrown
            // away (lambda/u_corr end the warmup already partway converged), and the estimate
            // reflects THIS step's actual active set, not a linearized approximation of it.
            // NOTE (see CONDENSED_SOLVER_REPORT.md): tracing this per-iteration on domino found the
            // observed ratio does NOT settle to a stable value even over 20 iterations -- it drifts
            // upward from ~0.3 to ~0.9 and transiently EXCEEDS 1 (residual growing) before this
            // warmup budget runs out, for a genuinely nonsmooth reason (the active contact set is
            // still resolving, not merely "not yet asymptotic"). This makes the observed-ratio
            // estimate unreliable regardless of warmup length, not just too-short-a-warmup -- kept
            // short (cheap) since a longer warmup measurably does NOT fix this.
            constexpr int kWarmupIters = 4;
            VectorXr u_prev = u_corr;
            real_t lastErr = -1, ratioEstimate = -1;
            for (int w = 0; w < kWarmupIters; ++w) {
                auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
                doSweep(lambda, u_corr);
                sc_sweep.~Scope();
                warmupItersUsed = w + 1;
                const real_t err = residualNorm(u_prev, u_corr);
                if (getenv("CARDILLO_DEBUG_WARMUP_TRACE")) std::cout << "  [warmup] w=" << w << " err=" << err << " ratio=" << (lastErr > 0 ? err / lastErr : -1) << std::endl;
                if (err <= (real_t)1) {
                    warmupConverged = true;  // converged during warmup alone -- nothing left to accelerate
                    break;
                }
                if (lastErr > (real_t)0 && std::isfinite(err) && std::isfinite(lastErr)) ratioEstimate = err / lastErr;
                lastErr = err;
                u_prev = u_corr;
            }
            if (!warmupConverged && ratioEstimate > (real_t)0 && std::isfinite(ratioEstimate)) rho = std::min(ratioEstimate, (real_t)0.995);
            if (m_cfg.debug_pj) {
                std::cout << "[Condensed] chebyshev adaptive rho: warmup_iters=" << warmupItersUsed << " ratio_estimate=" << ratioEstimate << " using_rho=" << rho
                          << (warmupConverged ? " (converged during warmup)" : "") << std::endl;
            }
        }

        VectorXr lambda_k = lambda, u_k = u_corr, lambda_y = lambda, u_y = u_corr;
        const double rho2 = (double)rho * (double)rho;
        double omega = 1.0;

        for (int iter = 0; !warmupConverged && iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            lambda = lambda_y;
            u_corr = u_y;
            doSweep(lambda, u_corr);
            sc_sweep.~Scope();

            const VectorXr lambda_k1 = lambda;
            const VectorXr u_k1 = u_corr;
            m_last_iters = warmupItersUsed + iter + 1;

            const real_t err = residualNorm(u_y, u_k1);
            if (m_cfg.debug_pj && iter % 1000 == 0) std::cout << "[Condensed] chebyshev iter " << iter << ", residual norm: " << err << std::endl;
            if (err <= (real_t)1) break;

            // If the RAW sweep output (lambda_k1/u_k1, before any extrapolation) is already
            // non-finite, the underlying sweep itself has diverged (e.g. condensed.sweep_mode=
            // jacobi at an alpha too high for its stability margin -- see condensed.rst's Nesterov
            // warning, the same pre-existing instability, not something this acceleration method
            // is expected to fix) -- there is nothing left to fall back to or extrapolate from, so
            // this must throw here, exactly like the Nesterov branch above's identical check.
            // Missing this check does NOT self-correct: falling back to lambda_y=lambda_k1 below
            // would just keep feeding the same NaN state into doSweep() every remaining iteration,
            // spinning silently for the rest of pj_max_iterations instead of failing fast.
            if (!lambda_k1.allFinite() || !u_k1.allFinite()) checkResidual(u_y, u_k1, iter);  // throws

            const double omega_next = (iter == 0) ? (1.0 / (1.0 - rho2 * 0.5)) : (1.0 / (1.0 - (rho2 * 0.25) * omega));

            if (!std::isfinite(omega_next) || omega_next <= 0.0 || !lambda_k1.allFinite() || !u_k1.allFinite()) {
                // The linearization this method relies on broke down (e.g. the active set changed
                // sharply due to projection) -- fall back to a single un-accelerated step rather
                // than risk amplifying a bad extrapolation. Matches PJ's own chebyshev_loop() guard.
                lambda_y = lambda_k1;
                u_y = u_k1;
                omega = 1.0;
            } else {
                // Guard the extrapolation itself, not just its aftermath -- same reasoning as the
                // identical guard in the Nesterov branch above: omega_next being finite/positive
                // bounds the *coefficient*, not the extrapolated magnitude.
                VectorXr trial_lambda_y = (real_t)omega_next * (lambda_k1 - lambda_k) + lambda_k;
                VectorXr trial_u_y = (real_t)omega_next * (u_k1 - u_k) + u_k;
                if (trial_lambda_y.allFinite() && trial_u_y.allFinite()) {
                    lambda_y = std::move(trial_lambda_y);
                    u_y = std::move(trial_u_y);
                    omega = omega_next;
                } else {
                    lambda_y = lambda_k1;
                    u_y = u_k1;
                    omega = 1.0;
                }
            }

            lambda_k = lambda_k1;
            u_k = u_k1;
        }
        // `lambda`/`u_corr` already hold the latest post-sweep values from inside the loop above --
        // no final reassignment needed, same fact the Nesterov branch above relies on. If the
        // adaptive-rho warmup itself already converged, the loop above never ran (its condition is
        // `!warmupConverged && ...`), so m_last_iters was never set by it -- set it here instead.
        if (warmupConverged) m_last_iters = warmupItersUsed;
    }

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) WarmstartProvider::storeImpulse(lambda.segment(nSprings + nDampers, nContacts), m_dyn, /*invertNormalSign=*/true);

    return u_free - u_corr;
}

}  // namespace cardillo::solver
