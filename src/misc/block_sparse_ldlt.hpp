#pragma once

#include <unordered_map>
#include <vector>
#include "types.hpp"

namespace cardillo::misc {

// Block-sparse LDLT (Cholesky-like) factorization of a symmetric matrix given as a graph of small
// dense blocks -- never assembles a dense or Eigen::SparseMatrix global matrix, matching the
// block-sparse philosophy of the condensed contact solver this was built for (see
// docs/chapters/solvers/condensed.rst's Schur-complement section). Nodes are 0..n-1, each with a
// fixed dof `dims[i]` (<=6 in every current use case).
//
// Two modes, selected at build() time:
//
// - `symmetric=true` (the default, and the only mode used before implicit-gyroscopic support):
//   off-diagonal coupling is given as an *undirected* edge list, symmetric by construction (edge
//   (i,j) implies block(j,i) = block(i,j)^T -- the caller supplies only one direction). Factors as
//   a true block LDLT (S = L D L^T), exploiting symmetry: one L list per pivot, half the Schur-
//   update work of the general case. This path is completely unchanged by the non-symmetric
//   addition below and must stay that way -- implicit gyroscopic forces are the only reason
//   anything here needs to be asymmetric at all, and most scenes never touch them.
// - `symmetric=false` (needed when implicit gyroscopic forces make a body's effective mass -- and
//   therefore any Sbb block touching it -- genuinely non-symmetric, block(i,j) != block(j,i)^T):
//   the edge list is *directed* -- the caller must supply both (i,j) and (j,i) explicitly whenever
//   the coupling is asymmetric (supplying only one direction leaves the other as zero). Factors as
//   a true block LU (S = L U, U not unit-triangular), storing both an L list (rows below) and a U
//   list (columns to the right) per pivot -- roughly double the Schur-update work of the symmetric
//   path, since nothing can be derived via transpose anymore. This is the whole reason the
//   symmetric path is kept separate rather than special-cased inside one general algorithm: paying
//   that cost on every scene, including the overwhelming majority that never enable implicit
//   gyroscopic forces, would be a pure regression for no benefit.
//
// Elimination order is never left to chance: factor() computes a greedy minimum-degree order
// internally (a lightweight symbolic pre-pass over adjacency only, no matrix blocks touched) --
// natural/creation order is a correct but NOT merely-suboptimal choice for a general graph: on a
// scene with real branching compliant networks (e.g. hangbridge's tripod/deck topology) it can
// cascade into a fill-in blowup that is not just slower but impractically so, confirmed by
// measurement (a from-scratch attempt at natural order hung on hangbridge). Minimum-degree is the
// standard, well-understood fix. For a simple path graph (every node degree <=2, e.g. every other
// current example scene's compliant network) minimum-degree recovers exactly the chain order with
// zero fill-in -- the same operation count as natural order there, so nothing is given up on the
// common case to fix the uncommon one. The same order-selection logic is reused unchanged for both
// modes -- it operates on adjacency only, never on whether the blocks happen to be symmetric.
class BlockSparseLDLT {
   public:
    // dims[i]: dof of node i. diagBlocks[i]: initial dims[i] x dims[i] diagonal block.
    // edgeNodes[e] = {i, j} with i != j, edgeBlocks[e] of shape dims[i] x dims[j]: the initial
    // (i,j) coupling.
    //   symmetric=true (default): edge list is undirected -- at most one entry per unordered pair;
    //     block(j,i) is derived as block(i,j)^T. diagBlocks must be symmetric (not required PD
    //     upfront -- invertSmallSpd() falls back gracefully per pivot).
    //   symmetric=false: edge list is directed -- (i,j) and (j,i) are independent entries; supply
    //     both explicitly whenever the true coupling is asymmetric (an omitted direction is taken
    //     as exactly zero, not "derive from the other direction"). diagBlocks need not be
    //     symmetric either.
    void build(std::vector<int> dims, std::vector<MatrixXXr> diagBlocks, const std::vector<std::array<int, 2>>& edgeNodes, const std::vector<MatrixXXr>& edgeBlocks, bool symmetric = true);

    // Factors in place using an internally-computed greedy minimum-degree elimination order. Call
    // once after build(), before solve(). Dispatches to the block-LDLT or block-LU algorithm
    // according to the `symmetric` flag passed to build().
    void factor();

    // Same factorization algorithm, but with a caller-supplied elimination order instead of the
    // internal minimum-degree heuristic. Exposed for testing correctness independent of ordering
    // choice (see the unit check this was validated with) -- production code should call factor().
    void factorWithOrder(std::vector<int> order);

    // Solves S * x = rhs, where S is the matrix passed to build() and rhs is stacked per-node in
    // node order 0..n-1 (sizes given by dims). Can be called repeatedly with different rhs after
    // one factor()/factorWithOrder() call.
    VectorXr solve(const VectorXr& rhs) const;

    int numNodes() const { return (int)m_dims.size(); }
    const std::vector<int>& dims() const { return m_dims; }
    int totalDim() const;

    // The elimination order actually used by the last factor()/factorWithOrder() call. Exposed so
    // a caller whose graph *structure* (dims + edges, not the numeric blocks) is stable across
    // repeated build()+factor() calls -- e.g. CondensedAssembler's bilateral graph, whose topology
    // doesn't change step-to-step even though Gii/complianceDiag do -- can cache it and pass it to
    // factorWithOrder() on a structural cache hit, skipping the O(n^2) minimum-degree symbolic pass
    // entirely. See CondensedAssembler::buildBilateralFactorization().
    const std::vector<int>& order() const { return m_order; }

   private:
    // Every node's dof is <=6 (see class comment), so every block below is a fixed-size Matrixr<6,6>
    // stack buffer accessed via `.topLeftCorner(dims[i], dims[j])`, exactly the Buf6-style pattern
    // used for the condensed solver's own per-block hot-path temporaries: this is the actual hot
    // loop (factor()'s Schur update touches every stored block once per pivot; solve()'s forward/
    // backward substitution touches every stored block once per call), unlike the one-time
    // ingestion in build(), so this is where fixing the storage size actually pays off.
    struct PivotFactor {
        int k{0};
        Matrixr<6, 6> Dinv{Matrixr<6, 6>::Zero()};                    // .topLeftCorner(dims[k], dims[k])
        std::vector<std::pair<int, Matrixr<6, 6>>> L;                 // (neighbor p, L_{p,k} = block(p,k)_at_elimination * Dinv), .topLeftCorner(dims[p], dims[k])
        // Non-symmetric mode only: U_{k,q} = block(k,q)_at_elimination (the raw current row-k
        // value, no Dinv applied -- Dinv is applied once, after subtracting U-contributions, in
        // the backward solve). Empty and unused in symmetric mode, where U is recovered from L via
        // transpose instead of stored separately.
        std::vector<std::pair<int, Matrixr<6, 6>>> U;                 // (neighbor q, U_{k,q}), .topLeftCorner(dims[k], dims[q])
    };

    // Symbolic-only simulation of the same fill-in process factorWithOrder() performs, tracking
    // just neighbor sets (not matrix blocks) -- cheap enough to run one full min-degree selection
    // pass upfront (O(n^2) node scans for the node counts seen in practice, hundreds to low
    // thousands) before ever touching the actual dense blocks. Shared by both modes -- fill-in
    // depends only on adjacency, not on whether the blocks are symmetric.
    std::vector<int> computeGreedyMinDegreeOrder() const;

    // The two factorization algorithms (block LDLT vs block LU), dispatched from factorWithOrder()
    // based on m_symmetric. Kept as fully separate functions rather than one function with branches
    // sprinkled through it, specifically so the symmetric path is trivially readable as "exactly
    // what it was before this mode existed" -- see the class comment.
    void factorSymmetric(const std::vector<int>& order);
    void factorNonSymmetric(const std::vector<int>& order);
    VectorXr solveSymmetric(const VectorXr& rhs) const;
    VectorXr solveNonSymmetric(const VectorXr& rhs) const;

    std::vector<int> m_dims;
    bool m_symmetric{true};
    std::vector<Matrixr<6, 6>> m_diag;                              // working diagonal blocks (.topLeftCorner(dims[i],dims[i])), mutated during factor()
    std::vector<std::unordered_map<int, Matrixr<6, 6>>> m_off;      // m_off[i][j]: block(i,j) (.topLeftCorner(dims[i],dims[j])). Symmetric mode: both directions stored, block(j,i) kept as its transpose. Non-symmetric mode: exactly what was supplied/computed for that direction, independent of m_off[j][i].
    std::vector<int> m_order;
    std::vector<PivotFactor> m_factors;  // filled by factor()/factorWithOrder(), in elimination order
};

}  // namespace cardillo::misc
