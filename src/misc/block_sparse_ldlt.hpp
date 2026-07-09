#pragma once

#include <unordered_map>
#include <vector>
#include "types.hpp"

namespace cardillo::misc {

// Block-sparse LDLT (Cholesky-like) factorization of a symmetric matrix given as a graph of small
// dense blocks -- never assembles a dense or Eigen::SparseMatrix global matrix, matching the
// block-sparse philosophy of the condensed contact solver this was built for (see
// docs/chapters/solvers/condensed.rst's Schur-complement section). Nodes are 0..n-1, each with a
// fixed dof `dims[i]` (<=6 in every current use case); off-diagonal coupling is given as a sparse
// edge list, symmetric by construction (edge (i,j) implies block(j,i) = block(i,j)^T).
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
// common case to fix the uncommon one.
class BlockSparseLDLT {
   public:
    // dims[i]: dof of node i. diagBlocks[i]: initial dims[i] x dims[i] diagonal block (must be
    // symmetric; not required PD upfront, invertSmallSpd() falls back gracefully per pivot).
    // edgeNodes[e] = {i, j} with i != j: the initial (i,j) coupling, edgeBlocks[e] of shape
    // dims[i] x dims[j]. At most one edge per unordered pair -- merge duplicates before calling.
    void build(std::vector<int> dims, std::vector<MatrixXXr> diagBlocks, const std::vector<std::array<int, 2>>& edgeNodes, const std::vector<MatrixXXr>& edgeBlocks);

    // Factors in place using an internally-computed greedy minimum-degree elimination order. Call
    // once after build(), before solve().
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
    };

    // Symbolic-only simulation of the same fill-in process factorWithOrder() performs, tracking
    // just neighbor sets (not matrix blocks) -- cheap enough to run one full min-degree selection
    // pass upfront (O(n^2) node scans for the node counts seen in practice, hundreds to low
    // thousands) before ever touching the actual dense blocks.
    std::vector<int> computeGreedyMinDegreeOrder() const;

    std::vector<int> m_dims;
    std::vector<Matrixr<6, 6>> m_diag;                              // working diagonal blocks (.topLeftCorner(dims[i],dims[i])), mutated during factor()
    std::vector<std::unordered_map<int, Matrixr<6, 6>>> m_off;      // m_off[i][j]: block(i,j) (.topLeftCorner(dims[i],dims[j])), both directions stored (block(j,i) kept as its transpose)
    std::vector<int> m_order;
    std::vector<PivotFactor> m_factors;  // filled by factor()/factorWithOrder(), in elimination order
};

}  // namespace cardillo::misc
