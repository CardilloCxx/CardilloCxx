#include "condensed_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>

#ifdef CARDILLO_HAVE_OPENMP
#include <omp.h>
#endif

#include "../../misc/contact_projection.hpp"
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
// views for whichever n <= 6 actually applies; unlike VectorXr it never touches the heap.
using Buf6 = Vector6r;

// Local update of a block from step = rho * xi (already relaxed/scaled), see the theory manual,
// eq:projection_update: bilateral blocks take the linear step lamOld - step, contact blocks the
// negative variant lambda = -proj_C(step - lamOld) (misc/contact_projection.hpp, shared with
// ProjectedGaussSeidel/ProjectedJacobi). Templated so it works on a VectorXr::segment() view or a
// Buf6::head() view. `o` holds the step on entry and the new lambda on exit.
template <typename Derived, typename DerivedOld>
void updateBlock(RowBlock::Kind kind, real_t mu, const Eigen::MatrixBase<DerivedOld>& lamOld, Eigen::MatrixBase<Derived>& o) {
    if (kind == RowBlock::Kind::ContactFrictionless)
        o[0] = misc::proxNormal(o[0] - lamOld[0]);
    else if (kind == RowBlock::Kind::ContactFrictional)
        misc::proxCoulomb(o - lamOld, mu, o);
    else
        o = lamOld - o;
}

// Writes this block's row velocity xi = xi^free + Ja*du_A + Jb*du_B + complianceDiag*lambda_block
// (theory manual, eq:block_residual; rhs = xi^free) into out.head(blk.dim). No heap allocation:
// `out` is caller-owned stack storage, written into via `.noalias()` rather than returned by value.
void blockResidual(const RowBlock& blk, const VectorXr& rhs, const VectorXr& du, const VectorXr& lambda, Buf6& out) {
    auto o = out.head(blk.dim);
    o = rhs.segment(blk.offset, blk.dim);
    if (!blk.bodyContributions.empty()) {
        for (const auto& contribution : blk.bodyContributions)
            o.noalias() += contribution.J * du.segment(contribution.offset, contribution.dof);
    } else {
        if (blk.aDof > 0) o.noalias() += blk.Ja * du.segment(blk.aOff, blk.aDof);
        if (blk.bDof > 0) o.noalias() += blk.Jb * du.segment(blk.bOff, blk.bDof);
    }
    o.noalias() += blk.complianceDiag.cwiseProduct(lambda.segment(blk.offset, blk.dim));
}

// nullptr for every body without an active implicit-gyroscopic override -- the overwhelming common
// case, so every call site below keeps its original plain-diagonal-MinvDiag expression unchanged in
// that branch. See CondensedTopology::gyroMinvBlocks.
const Matrix66r* gyroBlockFor(int bodyIndex, const std::unordered_map<int, Matrix66r>& gyroBlocks) {
    if (bodyIndex < 0 || gyroBlocks.empty()) return nullptr;
    auto it = gyroBlocks.find(bodyIndex);
    return it == gyroBlocks.end() ? nullptr : &it->second;
}

// Applies one block's impulse delta (dlambda.head(blk.dim)) to du in place: du += Minv *
// J^T * dlambda. `tmp` is caller-provided scratch for the intermediate J^T*dlambda product (also
// Buf6, also heap-free). `gyroBlocks` is topo.gyroMinvBlocks -- empty for every scene without an
// active implicit-gyroscopic body, in which case every branch below reduces to exactly the original
// expression (see gyroBlockFor()).
void scatterDelta(const RowBlock& blk, const VectorXr& MinvDiag, const std::unordered_map<int, Matrix66r>& gyroBlocks, const Buf6& dlambda, VectorXr& du, Buf6& tmp) {
    if (!blk.bodyContributions.empty()) {
        for (const auto& contribution : blk.bodyContributions) {
            tmp.head(contribution.dof).noalias() = contribution.J.transpose() * dlambda.head(blk.dim);
            if (const auto* gyro = gyroBlockFor(contribution.bodyIndex, gyroBlocks))
                du.segment(contribution.offset, contribution.dof).noalias() += gyro->topLeftCorner(contribution.dof, contribution.dof) * tmp.head(contribution.dof);
            else
                du.segment(contribution.offset, contribution.dof).noalias() += MinvDiag.segment(contribution.offset, contribution.dof).cwiseProduct(tmp.head(contribution.dof));
        }
        return;
    }
    if (blk.aDof > 0) {
        tmp.head(blk.aDof).noalias() = blk.Ja.transpose() * dlambda.head(blk.dim);
        if (const auto* gA = gyroBlockFor(blk.bodyIndexA, gyroBlocks))
            du.segment(blk.aOff, blk.aDof).noalias() += gA->topLeftCorner(blk.aDof, blk.aDof) * tmp.head(blk.aDof);
        else
            du.segment(blk.aOff, blk.aDof).noalias() += MinvDiag.segment(blk.aOff, blk.aDof).cwiseProduct(tmp.head(blk.aDof));
    }
    if (blk.bDof > 0) {
        tmp.head(blk.bDof).noalias() = blk.Jb.transpose() * dlambda.head(blk.dim);
        if (const auto* gB = gyroBlockFor(blk.bodyIndexB, gyroBlocks))
            du.segment(blk.bOff, blk.bDof).noalias() += gB->topLeftCorner(blk.bDof, blk.bDof) * tmp.head(blk.bDof);
        else
            du.segment(blk.bOff, blk.bDof).noalias() += MinvDiag.segment(blk.bOff, blk.bDof).cwiseProduct(tmp.head(blk.bDof));
    }
}

// condensed.true_schur=true: exactly zeroes every bilateral (spring+damper) row's residual in one
// shot, given the CURRENT contact contribution to du, via the precomputed block-sparse LDLT of
// Sbb (see CondensedAssembler::buildBilateralFactorization()). This replaces those rows'
// participation in the Gauss-Seidel/Jacobi/colored sweep entirely -- bilateral rows are pure
// equality constraints (updateBlock() takes the plain linear step for them), so unlike the projected contact
// rows there is no iteration needed once their joint linear system is solved exactly. A no-op when
// there are no bilateral blocks (e.g. domino), by construction.
void exactBilateralStep(const CondensedTopology& topo, const misc::BlockSparseLDLT& ldlt, const VectorXr& rhs, const VectorXr& MinvDiag, VectorXr& lambda, VectorXr& du) {
    const int n = topo.numBilateralBlocks;
    if (n == 0) return;

    VectorXr xi(ldlt.totalDim());
    Buf6 rb;
    for (int i = 0; i < n; ++i) {
        const auto& blk = topo.blocks[(size_t)i];
        blockResidual(blk, rhs, du, lambda, rb);
        xi.segment(blk.offset, blk.dim) = rb.head(blk.dim);
    }

    // pi_b <- pi_b - Sbb^{-1} xi_b (theory manual, eq:exact_bilateral_step)
    const VectorXr delta = -ldlt.solve(xi);

    Buf6 d, tmp;
    for (int i = 0; i < n; ++i) {
        const auto& blk = topo.blocks[(size_t)i];
        d.head(blk.dim) = delta.segment(blk.offset, blk.dim);
        scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, d, du, tmp);
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
// branch and the projection in updateBlock(), making this update a pure linear relaxation step
// (o = lamOld - relaxation*alpha_i*GiiInv*xi, nothing else). Used only by estimateSpectralRadius()
// below, which needs the *linear* part of the sweep operator in isolation -- Newton/projection are
// exactly the nonlinearities that make a rigorous spectral-radius argument inapplicable in the
// first place (see that function's own comment).
void computeBlockUpdate(const RowBlock& blk, const Buf6& xi, real_t relaxation, real_t alpha, const Buf6& lamOld, bool useNewton, const NewtonACParams& newtonParams, Buf6& out,
                         bool linearOnly = false) {
    if (useNewton && !linearOnly && blk.kind == RowBlock::Kind::ContactFrictional) {
        Vector3r lam3 = lamOld.head<3>();
        const Vector3r xi3 = xi.head<3>();
        if (solveContactBlockNewtonAC(blk.Gii, xi3, blk.mu, lam3, newtonParams)) {
            out.head<3>() = lam3;
            return;
        }
    }

    // step = sigma * alpha_i * P_i * xi
    auto o = out.head(blk.dim);
    o.noalias() = relaxation * (blk.GiiInv * xi.head(blk.dim));
    if (isContactKind(blk.kind)) o *= alpha;
    if (linearOnly)
        o = lamOld.head(blk.dim) - o;  // updateBlock() without projection
    else
        updateBlock(blk.kind, blk.mu, lamOld.head(blk.dim), o);
}

// Sequential block-Gauss-Seidel sweep: propagates du immediately after each block, matching
// ProjectedGaussSeidelSolver's own loop structure exactly, just expressed over RowBlocks.
// `startBlock` skips the bilateral (spring+damper) prefix when condensed.true_schur handles those
// exactly instead (0 -- i.e. every block -- otherwise); contact blocks are always a contiguous
// suffix of topo.blocks, so a single start index is enough, no separate index list needed.
void gaussSeidelSweep(const CondensedTopology& topo, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton, const NewtonACParams& newtonParams,
                       VectorXr& lambda, VectorXr& du, int startBlock = 0, bool linearOnly = false) {
    Buf6 r, lamOld, lamNew, delta, tmp;
    for (size_t bi = (size_t)startBlock; bi < topo.blocks.size(); ++bi) {
        const auto& blk = topo.blocks[bi];
        blockResidual(blk, rhs, du, lambda, r);
        lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
        computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew, linearOnly);

        delta.head(blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
        scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, delta, du, tmp);
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
// writes of du touch disjoint index ranges.
void coloredSweep(const CondensedTopology& topo, const misc::Coloring& coloring, const VectorXr& rhs, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, bool useNewton,
                   const NewtonACParams& newtonParams, VectorXr& lambda, VectorXr& du, bool linearOnly = false) {
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
            blockResidual(blk, rhs, du, lambda, r);
            lamOld.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
            computeBlockUpdate(blk, r, relaxation, alpha, lamOld, useNewton, newtonParams, lamNew, linearOnly);

            // Safe: no two blocks in `cls` share a body (that's the coloring invariant), so these
            // writes land in disjoint slices of du/lambda across threads -- no race, no atomics.
            delta.head(blk.dim) = lamNew.head(blk.dim) - lamOld.head(blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, delta, du, tmp);
            lambda.segment(blk.offset, blk.dim) = lamNew.head(blk.dim);
        }
        // implicit barrier at the end of `omp for` here, before the next color starts
    }
}

// Estimates the spectral radius of the LINEARIZED sweep operator (i.e. with projection/Newton
// disabled -- exactly the nonlinearities that make a rigorous spectral-radius argument
// inapplicable in the first place, same caveat as ProjectedJacobiSolver's own estimator), for
// Chebyshev semi-iterative acceleration (condensed.pj_chebyshev... see pj.chebyshev). Mirrors
// ProjectedJacobiSolver's own estimateSpectralRadius(), adapted to condensed's matrix-free
// block-sweep representation instead of a global sparse system: rather than needing an explicit
// linear operator/solveS(), this runs the SAME sweep function the real iteration uses
// (gaussSeidelSweep/coloredSweep, plus exactBilateralStep under true_schur) with
// projection/Newton disabled (linearOnly=true) and a ZERO rhs. With no projection and no rhs, one
// sweep step is an EXACT linear map G(d) = A*d (no affine offset -- the rhs-dependent constant
// term vanishes), so plain power iteration on a probe vector converges to |A|'s dominant
// eigenvalue magnitude, exactly analogous to PJ's own power iteration on its (implicit,
// matrix-based) linearized operator. Uses a deterministic all-ones probe vector (matching PJ's own
// choice) rather than a randomized one, so this doesn't introduce a new nondeterminism source into
// the solver.
real_t estimateSpectralRadius(const CondensedTopology& topo, const VectorXr& MinvDiag, real_t relaxation, real_t alpha, const misc::Coloring* coloring, int startBlock, bool useTrueSchur,
                              const misc::BlockSparseLDLT& bilateralLdlt, int nV) {
    constexpr int kPowerIterations = 8;
    const int numLambda = topo.numLambda;
    if (numLambda == 0) return (real_t)0;

    const VectorXr zeroRhs = VectorXr::Zero(numLambda);
    const NewtonACParams unusedNewtonParams{};  // never consulted: linearOnly=true skips the Newton branch entirely

    VectorXr d_lambda = VectorXr::Ones(numLambda).normalized();
    VectorXr d_du = VectorXr::Zero(nV);
    {
        Buf6 lam0, tmp0;
        for (const auto& blk : topo.blocks) {
            lam0.head(blk.dim) = d_lambda.segment(blk.offset, blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, lam0, d_du, tmp0);
        }
    }

    real_t rho = 0;
    for (int it = 0; it < kPowerIterations; ++it) {
        if (useTrueSchur) exactBilateralStep(topo, bilateralLdlt, zeroRhs, MinvDiag, d_lambda, d_du);
        if (coloring) {
            coloredSweep(topo, *coloring, zeroRhs, MinvDiag, relaxation, alpha, false, unusedNewtonParams, d_lambda, d_du, /*linearOnly=*/true);
        } else {
            gaussSeidelSweep(topo, zeroRhs, MinvDiag, relaxation, alpha, false, unusedNewtonParams, d_lambda, d_du, startBlock, /*linearOnly=*/true);
        }

        const real_t n = d_lambda.norm();
        // When true_schur solves the whole probe exactly (no contacts), n is round-off of an O(1)
        // computation, not an eigenvalue; renormalizing it would amplify noise to a spurious rho ~ 1.
        if (!(n > (real_t)1e-8)) return (real_t)0;
        rho = n;
        d_lambda /= n;
        d_du /= n;
    }
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

    // Warmstart -- identical boundary behavior to PGS (contact multipliers are the canonical
    // percussions, see warmstart.hpp).
    VectorXr lambda = VectorXr::Zero(topo.numLambda);
    {
        VectorXr Lambda_g = m_dyn.Lambda_g();
        if ((int)Lambda_g.size() == nSprings) lambda.head(nSprings) = Lambda_g;
        VectorXr Lambda_gamma = m_dyn.Lambda_gamma();
        if ((int)Lambda_gamma.size() == nDampers) lambda.segment(nSprings, nDampers) = Lambda_gamma;
    }
    if (m_cfg.pj_warmstart && nContacts > 0) {
        VectorXr l_contact = VectorXr::Zero(nContacts);
        WarmstartProvider::applyWarmstart(l_contact, m_dyn);
        lambda.segment(nSprings + nDampers, nContacts) = l_contact;
    }

    VectorXr du = VectorXr::Zero(nV);
    {
        Buf6 lam0, tmp0;
        for (const auto& blk : topo.blocks) {
            lam0.head(blk.dim) = lambda.segment(blk.offset, blk.dim);
            scatterDelta(blk, MinvDiag, topo.gyroMinvBlocks, lam0, du, tmp0);
        }
    }

#ifdef CARDILLO_HAVE_OPENMP
    if (m_cfg.condensed_num_threads > 0) omp_set_num_threads(m_cfg.condensed_num_threads);
#endif

    sc_setup.stop();
    auto sc_solve = m_dyn.timings()->scope(misc::TimingManager::TimerId::Condensed);

    const real_t alpha = m_cfg.pj_alpha;
    const real_t relaxation = m_cfg.pj_relaxation;

    VectorXr res_v0(nV), res_v1(nV), res_diff(nV), res_q(nV);
    auto residualNorm = [&](const VectorXr& u_prev, const VectorXr& u_cur) -> real_t {
        res_diff.noalias() = u_cur - u_prev;
        res_v0.noalias() = u_free + u_prev;
        res_v1.noalias() = u_free + u_cur;
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

    // condensed.true_schur eliminates the bilateral (spring+damper) blocks exactly in every sweep;
    // the projected sweep then starts at the first contact block.
    const bool useTrueSchur = m_cfg.condensed_true_schur;
    const int startBlock = useTrueSchur ? topo.numBilateralBlocks : 0;
    const bool colored = m_cfg.condensed_sweep_mode == "colored";

    misc::BlockSparseLDLT bilateralLdlt;
    if (useTrueSchur && topo.numBilateralBlocks > 0) {
        auto sc_schur = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);
        bilateralLdlt = m_assembler.buildBilateralFactorization(topo);
    }

    misc::Coloring coloring;
    if (colored) {
        auto sc_color = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedColoring);
        const auto adjacency = buildBlockAdjacency(topo, startBlock);
        coloring = misc::colorGreedyWelshPowell((int)topo.blocks.size() - startBlock, adjacency);
        for (auto& cls : coloring.colorClasses)
            for (auto& idx : cls) idx += startBlock;
    }

    const bool useNewton = (m_cfg.condensed_local_solve == "newton");
    const NewtonACParams newtonParams{m_cfg.condensed_newton_max_iters, m_cfg.condensed_newton_tol,
                                       (m_cfg.condensed_newton_rho_strategy == "full") ? NewtonRhoStrategy::FullSpectral : NewtonRhoStrategy::Split};

    // pj.chebyshev and pj.nesterov are shared with PJ/PGS, with the same precedence (Nesterov wins).
    const bool useChebyshev = m_cfg.pj_chebyshev && !m_cfg.pj_nesterov;
    real_t chebyshevRho = 0;
    if (useChebyshev) {
        // Chebyshev's spectral-radius estimate is only rigorously meaningful when the linearized
        // sweep operator has a REAL eigenvalue spectrum -- guaranteed (via similarity to a
        // symmetric matrix) when Minv is symmetric, i.e. no active implicit-gyroscopic override
        // anywhere in the system (see docs/chapters/solvers/condensed.rst's Chebyshev section).
        // A non-symmetric Minv can give the operator complex eigenvalues, which the classical
        // 3-term Chebyshev recurrence isn't designed for -- not a guaranteed failure (the same
        // "heuristic on the full nonlinear problem" caveat already applies regardless), but worth
        // surfacing rather than silently trusting the estimate.
        if (m_cfg.debug_pj && !topo.gyroMinvBlocks.empty()) {
            std::cout << "[Condensed] Warning: pj.chebyshev's spectral-radius estimate assumes a real "
                         "eigenvalue spectrum (via symmetric Minv); an implicit-gyroscopic body is active, "
                         "making Minv non-symmetric -- the estimate may not be reliable here."
                      << std::endl;
        }
        auto sc_cheb = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSetup);
        const real_t rawRho = estimateSpectralRadius(topo, MinvDiag, relaxation, alpha, colored ? &coloring : nullptr, startBlock, useTrueSchur, bilateralLdlt, nV);
        // Clamp strictly below 1: the omega recurrence divides by (1 - rho^2/4*omega).
        chebyshevRho = std::min(rawRho, (real_t)0.995);
        if (m_cfg.debug_pj) std::cout << "[Condensed] chebyshev spectral radius estimate: " << chebyshevRho << std::endl;
    }

    // One sweep over (lambda, du); Nesterov and Chebyshev extrapolate between calls. Under
    // condensed.true_schur every call first solves the bilateral rows exactly for the current
    // contact contribution to du, then sweeps the contact blocks only.
    auto doSweep = [&](VectorXr& lam, VectorXr& ucorr) {
        if (useTrueSchur) exactBilateralStep(topo, bilateralLdlt, rhs, MinvDiag, lam, ucorr);
        if (colored) {
            coloredSweep(topo, coloring, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lam, ucorr);
        } else {
            gaussSeidelSweep(topo, rhs, MinvDiag, relaxation, alpha, useNewton, newtonParams, lam, ucorr, startBlock);
        }
    };

    if (!m_cfg.pj_nesterov && !useChebyshev) {
        VectorXr u_prev = du;
        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            doSweep(lambda, du);
            sc_sweep.stop();

            m_last_iters = iter + 1;
            if (checkResidual(u_prev, du, iter) < (real_t)1) break;
            u_prev = du;
        }
    } else if (m_cfg.pj_nesterov) {
        // Nesterov (FISTA-style) acceleration, ported verbatim from
        // ProjectedGaussSeidelSolver::solve()'s Nesterov branch: momentum-extrapolate between
        // sweeps, with adaptive restart when the residual grows, the momentum coefficient gets too
        // large, or the update direction reverses. Reuses the same shared config keys PJ/PGS
        // already use (pj.nesterov, pj.nesterov_beta_threshold, pj.nesterov_restart_limit) rather
        // than inventing condensed-specific ones.
        VectorXr lambda_k = lambda, u_k = du, lambda_y = lambda, u_y = du;
        double thk = 1.0;
        real_t err_prev = std::numeric_limits<real_t>::infinity();
        int restarts = 0;
        bool momentum_disabled = false;

        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            lambda = lambda_y;
            du = u_y;
            doSweep(lambda, du);
            sc_sweep.stop();

            const VectorXr lambda_k1 = lambda;
            const VectorXr u_k1 = du;
            m_last_iters = iter + 1;

            // A non-finite err is a restart trigger below, not an immediate failure.
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
        // `lambda`/`du` already hold the latest post-sweep values from inside the loop above
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
        const real_t rho = chebyshevRho;
        VectorXr lambda_k = lambda, u_k = du, lambda_y = lambda, u_y = du;
        const double rho2 = (double)rho * (double)rho;
        double omega = 1.0;

        for (int iter = 0; iter < m_cfg.pj_max_iterations; ++iter) {
            auto sc_sweep = m_dyn.timings()->scope(misc::TimingManager::TimerId::CondensedSweep);
            lambda = lambda_y;
            du = u_y;
            doSweep(lambda, du);
            sc_sweep.stop();

            const VectorXr lambda_k1 = lambda;
            const VectorXr u_k1 = du;
            m_last_iters = iter + 1;

            const real_t err = residualNorm(u_y, u_k1);
            if (m_cfg.debug_pj && iter % 1000 == 0) std::cout << "[Condensed] chebyshev iter " << iter << ", residual norm: " << err << std::endl;
            if (err <= (real_t)1) break;

            // If the RAW sweep output (lambda_k1/u_k1, before any extrapolation) is already
            // non-finite, the sweep itself has diverged and there is nothing to fall back to.
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
    }

    // The contact sweep of the last iteration changed du after its exact bilateral step. One more
    // step makes the returned velocities satisfy the bilateral rows exactly for the final contact
    // percussions, like those of projected Jacobi.
    if (useTrueSchur && nContacts > 0) exactBilateralStep(topo, bilateralLdlt, rhs, MinvDiag, lambda, du);

    if (nSprings > 0) m_dyn.setLambda_g(lambda.head(nSprings));
    if (nDampers > 0) m_dyn.setLambda_gamma(lambda.segment(nSprings, nDampers));
    if (nContacts > 0) WarmstartProvider::storeImpulse(lambda.segment(nSprings + nDampers, nContacts), m_dyn);

    return u_free + du;
}

}  // namespace cardillo::solver
