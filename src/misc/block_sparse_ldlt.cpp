#include "block_sparse_ldlt.hpp"

#include <limits>
#include <set>

#include "block_diagonal.hpp"  // invertSmallSpd, invertSmallGeneral

namespace cardillo::misc {

void BlockSparseLDLT::build(std::vector<int> dims, std::vector<MatrixXXr> diagBlocks, const std::vector<std::array<int, 2>>& edgeNodes, const std::vector<MatrixXXr>& edgeBlocks, bool symmetric) {
    m_dims = std::move(dims);
    m_symmetric = symmetric;
    const int n = (int)m_dims.size();

    m_diag.assign(n, Matrixr<6, 6>::Zero());
    for (int i = 0; i < n; ++i) m_diag[(size_t)i].topLeftCorner(m_dims[(size_t)i], m_dims[(size_t)i]) = diagBlocks[(size_t)i];

    m_off.assign(n, {});
    m_factors.clear();

    if (m_symmetric) {
        // Undirected: caller supplies one direction, the other is derived as its transpose.
        for (size_t e = 0; e < edgeNodes.size(); ++e) {
            const int i = edgeNodes[e][0];
            const int j = edgeNodes[e][1];
            Matrixr<6, 6> Bij = Matrixr<6, 6>::Zero();
            Bij.topLeftCorner(m_dims[(size_t)i], m_dims[(size_t)j]) = edgeBlocks[e];
            Matrixr<6, 6> Bji = Matrixr<6, 6>::Zero();
            Bji.topLeftCorner(m_dims[(size_t)j], m_dims[(size_t)i]) = edgeBlocks[e].transpose();
            m_off[(size_t)i][j] = Bij;
            m_off[(size_t)j][i] = Bji;
        }
    } else {
        // Directed: each entry is exactly one direction, stored as given -- the caller must supply
        // both (i,j) and (j,i) explicitly if both directions are nonzero (an omitted direction is
        // simply absent from m_off, i.e. treated as zero).
        for (size_t e = 0; e < edgeNodes.size(); ++e) {
            const int i = edgeNodes[e][0];
            const int j = edgeNodes[e][1];
            Matrixr<6, 6> Bij = Matrixr<6, 6>::Zero();
            Bij.topLeftCorner(m_dims[(size_t)i], m_dims[(size_t)j]) = edgeBlocks[e];
            m_off[(size_t)i][j] = Bij;
        }
    }
}

int BlockSparseLDLT::totalDim() const {
    int total = 0;
    for (int d : m_dims) total += d;
    return total;
}

std::vector<int> BlockSparseLDLT::computeGreedyMinDegreeOrder() const {
    // Adjacency for ordering purposes is symmetrized regardless of mode: in the non-symmetric case
    // the caller-supplied directed edges are still structurally symmetric in practice (asymmetry
    // is in the *values*, e.g. a gyroscopic-corrected mass block, not in *which* bodies/rows are
    // coupled -- see BlockSparseLDLT's class comment), so treating (i,j) or (j,i) present as "i and
    // j are adjacent" for fill-in purposes is correct either way.
    const int n = (int)m_dims.size();
    std::vector<std::set<int>> adj(n);
    for (int i = 0; i < n; ++i) {
        for (const auto& kv : m_off[(size_t)i]) {
            adj[(size_t)i].insert(kv.first);
            adj[(size_t)kv.first].insert(i);
        }
    }

    std::vector<bool> eliminated(n, false);
    std::vector<int> order;
    order.reserve(n);

    for (int step = 0; step < n; ++step) {
        int best = -1;
        size_t bestDeg = std::numeric_limits<size_t>::max();
        for (int i = 0; i < n; ++i) {
            if (eliminated[i]) continue;
            const size_t deg = adj[(size_t)i].size();
            if (deg < bestDeg) {
                bestDeg = deg;
                best = i;
            }
        }
        order.push_back(best);
        eliminated[best] = true;

        // Fill-in: every pair of `best`'s still-active neighbors becomes connected (mirrors the
        // real Schur update in factorWithOrder(), tracked here as adjacency only).
        std::vector<int> active;
        for (int nb : adj[(size_t)best]) {
            if (!eliminated[nb]) active.push_back(nb);
        }
        for (size_t a = 0; a < active.size(); ++a) {
            for (size_t b = a + 1; b < active.size(); ++b) {
                adj[(size_t)active[a]].insert(active[b]);
                adj[(size_t)active[b]].insert(active[a]);
            }
        }
        for (int nb : adj[(size_t)best]) adj[(size_t)nb].erase(best);
    }
    return order;
}

void BlockSparseLDLT::factor() { factorWithOrder(computeGreedyMinDegreeOrder()); }

void BlockSparseLDLT::factorWithOrder(std::vector<int> order) {
    m_order = std::move(order);
    if (m_symmetric) {
        factorSymmetric(m_order);
    } else {
        factorNonSymmetric(m_order);
    }
}

// Unchanged from before the non-symmetric mode existed -- exploits block(k,q) = block(q,k)^T, so
// only one L list per pivot is needed (U is recovered from L via transpose in solveSymmetric()).
void BlockSparseLDLT::factorSymmetric(const std::vector<int>& order) {
    m_factors.clear();
    m_factors.reserve(m_dims.size());
    std::vector<bool> eliminated(m_dims.size(), false);

    for (int k : order) {
        const int dimK = m_dims[(size_t)k];
        PivotFactor pf;
        pf.k = k;
        pf.Dinv.topLeftCorner(dimK, dimK) = invertSmallSpd(MatrixXXr(m_diag[(size_t)k].topLeftCorner(dimK, dimK)));

        // Still-active neighbors of k (copied out before mutating anything -- the update loop
        // below writes into m_off[p]/m_off[q] for p,q in this list, never into m_off[k] itself).
        std::vector<int> nbrs;
        nbrs.reserve(m_off[(size_t)k].size());
        for (const auto& kv : m_off[(size_t)k]) {
            if (!eliminated[(size_t)kv.first]) nbrs.push_back(kv.first);
        }

        pf.L.reserve(nbrs.size());
        for (int p : nbrs) {
            const int dimP = m_dims[(size_t)p];
            // Uninitialized, not Zero(): only ever read back via the same .topLeftCorner(dimP,dimK)
            // written here (dims[] are fixed for the object's lifetime), so the padding is never
            // read -- skip the wasted zero-fill of the untouched region.
            Matrixr<6, 6> Lpk;
            Lpk.topLeftCorner(dimP, dimK).noalias() = m_off[(size_t)p].at(k).topLeftCorner(dimP, dimK) * pf.Dinv.topLeftCorner(dimK, dimK);
            pf.L.emplace_back(p, std::move(Lpk));
        }

        // Right-looking Schur update over every pair (p, q) of still-active neighbors, including
        // p==q (the diagonal update) -- this is the block generalization of one elimination step
        // of a scalar sparse Cholesky/Thomas algorithm. Creates fill-in edges for p!=q pairs not
        // already connected.
        for (size_t ii = 0; ii < nbrs.size(); ++ii) {
            const int p = nbrs[ii];
            const int dimP = m_dims[(size_t)p];
            const auto Lp = pf.L[ii].second.topLeftCorner(dimP, dimK);  // dims[p] x dims[k]
            for (size_t jj = ii; jj < nbrs.size(); ++jj) {
                const int q = nbrs[jj];
                const int dimQ = m_dims[(size_t)q];
                const auto Bkq = m_off[(size_t)k].at(q).topLeftCorner(dimK, dimQ);  // dims[k] x dims[q]
                // Write straight into the destination (dims[p] x dims[q]) -- no padded 6x6
                // intermediate: Lp/Bkq are already exact-sized views, so `Lp * Bkq` evaluates
                // directly into the (correctly dims[p] x dims[q]-sized) `.topLeftCorner()`
                // destination via `.noalias() -=`, same operation count as the pre-fixed-size code.
                if (p == q) {
                    auto Dp = m_diag[(size_t)p].topLeftCorner(dimP, dimP);
                    Dp.noalias() -= Lp * Bkq;
                    Dp = (real_t)0.5 * (Dp + Dp.transpose());  // keep numerically symmetric
                } else {
                    auto it = m_off[(size_t)p].find(q);
                    if (it == m_off[(size_t)p].end()) {
                        // Fresh fill-in edge: no prior value to subtract from, so assign the
                        // negation directly rather than zero-filling first then subtracting.
                        Matrixr<6, 6> fresh;
                        fresh.topLeftCorner(dimP, dimQ).noalias() = -(Lp * Bkq);
                        it = m_off[(size_t)p].emplace(q, std::move(fresh)).first;
                    } else {
                        it->second.topLeftCorner(dimP, dimQ).noalias() -= Lp * Bkq;
                    }
                    Matrixr<6, 6> transposed;
                    transposed.topLeftCorner(dimQ, dimP) = it->second.topLeftCorner(dimP, dimQ).transpose();
                    m_off[(size_t)q][p] = std::move(transposed);
                }
            }
        }

        eliminated[(size_t)k] = true;
        m_factors.push_back(std::move(pf));
    }
}

// General block LU (S = L U, U not unit-triangular): block(k,q) and block(q,k) are independent, so
// both an L list (rows below, for forward substitution) and a U list (columns to the right, for
// backward substitution) are stored per pivot -- nothing can be derived from the other via
// transpose anymore. Roughly double factorSymmetric()'s Schur-update work: every unordered pair
// (p,q) of active neighbors needs *two* independent block updates (p,q) and (q,p), not one plus its
// transpose. This is the whole point of keeping the two paths separate rather than merging them --
// see the class comment.
void BlockSparseLDLT::factorNonSymmetric(const std::vector<int>& order) {
    m_factors.clear();
    m_factors.reserve(m_dims.size());
    std::vector<bool> eliminated(m_dims.size(), false);

    for (int k : order) {
        const int dimK = m_dims[(size_t)k];
        PivotFactor pf;
        pf.k = k;
        pf.Dinv.topLeftCorner(dimK, dimK) = invertSmallGeneral(MatrixXXr(m_diag[(size_t)k].topLeftCorner(dimK, dimK)));

        // Active neighbors from *either* direction: p is a neighbor if block(p,k) and/or block(k,p)
        // exists. Structurally the two directions have the same support in every current use case
        // (see the class comment), but iterate the union defensively rather than assume it.
        std::set<int> nbrSet;
        for (const auto& kv : m_off[(size_t)k])
            if (!eliminated[(size_t)kv.first]) nbrSet.insert(kv.first);
        for (int cand = 0; cand < (int)m_dims.size(); ++cand) {
            if (eliminated[(size_t)cand] || cand == k) continue;
            auto it = m_off[(size_t)cand].find(k);
            if (it != m_off[(size_t)cand].end()) nbrSet.insert(cand);
        }
        std::vector<int> nbrs(nbrSet.begin(), nbrSet.end());

        auto blockOrZero = [&](int a, int b, int dimA, int dimB) -> Matrixr<6, 6> {
            Matrixr<6, 6> out = Matrixr<6, 6>::Zero();
            auto it = m_off[(size_t)a].find(b);
            if (it != m_off[(size_t)a].end()) out.topLeftCorner(dimA, dimB) = it->second.topLeftCorner(dimA, dimB);
            return out;
        };

        // L_{p,k} = block(p,k) * Dinv_k (rows below), U_{k,q} = block(k,q) (columns to the right,
        // stored raw -- Dinv is applied once, after the backward-substitution subtraction, not
        // here; see solveNonSymmetric()).
        pf.L.reserve(nbrs.size());
        pf.U.reserve(nbrs.size());
        for (int p : nbrs) {
            const int dimP = m_dims[(size_t)p];
            Matrixr<6, 6> Apk = blockOrZero(p, k, dimP, dimK);
            Matrixr<6, 6> Lpk;
            Lpk.topLeftCorner(dimP, dimK).noalias() = Apk.topLeftCorner(dimP, dimK) * pf.Dinv.topLeftCorner(dimK, dimK);
            pf.L.emplace_back(p, std::move(Lpk));

            Matrixr<6, 6> Ukq = blockOrZero(k, p, dimK, dimP);
            pf.U.emplace_back(p, std::move(Ukq));
        }

        // Schur update over every ORDERED pair (p, q) of active neighbors (including p==q, the
        // diagonal update): block(p,q) -= L_{p,k} * U_{k,q}. Unlike the symmetric path, (p,q) and
        // (q,p) are independent and both must be computed directly -- neither is the other's
        // transpose in general.
        for (size_t ii = 0; ii < nbrs.size(); ++ii) {
            const int p = nbrs[ii];
            const int dimP = m_dims[(size_t)p];
            const auto Lp = pf.L[ii].second.topLeftCorner(dimP, dimK);
            for (size_t jj = 0; jj < nbrs.size(); ++jj) {
                const int q = nbrs[jj];
                const int dimQ = m_dims[(size_t)q];
                const auto Ukq = pf.U[jj].second.topLeftCorner(dimK, dimQ);

                if (p == q) {
                    m_diag[(size_t)p].topLeftCorner(dimP, dimP).noalias() -= Lp * Ukq;
                    continue;
                }
                auto it = m_off[(size_t)p].find(q);
                if (it == m_off[(size_t)p].end()) {
                    Matrixr<6, 6> fresh;
                    fresh.topLeftCorner(dimP, dimQ).noalias() = -(Lp * Ukq);
                    m_off[(size_t)p].emplace(q, std::move(fresh));
                } else {
                    it->second.topLeftCorner(dimP, dimQ).noalias() -= Lp * Ukq;
                }
            }
        }

        eliminated[(size_t)k] = true;
        m_factors.push_back(std::move(pf));
    }
}

VectorXr BlockSparseLDLT::solve(const VectorXr& rhs) const { return m_symmetric ? solveSymmetric(rhs) : solveNonSymmetric(rhs); }

VectorXr BlockSparseLDLT::solveSymmetric(const VectorXr& rhs) const {
    const int n = (int)m_dims.size();
    std::vector<int> offset(n + 1, 0);
    for (int i = 0; i < n; ++i) offset[(size_t)i + 1] = offset[(size_t)i] + m_dims[(size_t)i];

    std::vector<Vectorr<6>> y(n, Vectorr<6>::Zero());
    for (int i = 0; i < n; ++i) y[(size_t)i].head(m_dims[(size_t)i]) = rhs.segment(offset[(size_t)i], m_dims[(size_t)i]);

    // Forward substitution (L y = rhs, L unit lower-triangular in elimination order): visiting
    // pivots in elimination order guarantees y[pf.k] already holds every correction from earlier
    // pivots before it propagates its own effect onto its (later-eliminated) neighbors.
    for (const auto& pf : m_factors) {
        const int dimK = m_dims[(size_t)pf.k];
        for (const auto& pr : pf.L) {
            const int dimP = m_dims[(size_t)pr.first];
            y[(size_t)pr.first].head(dimP).noalias() -= pr.second.topLeftCorner(dimP, dimK) * y[(size_t)pf.k].head(dimK);
        }
    }

    std::vector<Vectorr<6>> x = y;
    for (const auto& pf : m_factors) {
        const int dimK = m_dims[(size_t)pf.k];
        x[(size_t)pf.k].head(dimK) = pf.Dinv.topLeftCorner(dimK, dimK) * y[(size_t)pf.k].head(dimK);
    }

    // Backward substitution (L^T x = D^-1 y), reverse elimination order.
    for (auto it = m_factors.rbegin(); it != m_factors.rend(); ++it) {
        const int dimK = m_dims[(size_t)it->k];
        for (const auto& pr : it->L) {
            const int dimP = m_dims[(size_t)pr.first];
            x[(size_t)it->k].head(dimK).noalias() -= pr.second.topLeftCorner(dimP, dimK).transpose() * x[(size_t)pr.first].head(dimP);
        }
    }

    VectorXr out(offset[(size_t)n]);
    for (int i = 0; i < n; ++i) out.segment(offset[(size_t)i], m_dims[(size_t)i]) = x[(size_t)i].head(m_dims[(size_t)i]);
    return out;
}

// S = L U (U not unit-triangular): forward substitution for L is identical in structure to the
// symmetric case (L is unit-lower-triangular either way). Backward substitution differs: solves
// U x = y directly using the stored raw U_{k,q} blocks, applying Dinv once *after* subtracting the
// U-contributions (rather than applying it up front and using L^T, which only works when U = D L^T,
// i.e. only in the symmetric case).
VectorXr BlockSparseLDLT::solveNonSymmetric(const VectorXr& rhs) const {
    const int n = (int)m_dims.size();
    std::vector<int> offset(n + 1, 0);
    for (int i = 0; i < n; ++i) offset[(size_t)i + 1] = offset[(size_t)i] + m_dims[(size_t)i];

    std::vector<Vectorr<6>> y(n, Vectorr<6>::Zero());
    for (int i = 0; i < n; ++i) y[(size_t)i].head(m_dims[(size_t)i]) = rhs.segment(offset[(size_t)i], m_dims[(size_t)i]);

    // Forward substitution (L y = rhs) -- same structure as the symmetric case, L is unit
    // lower-triangular in elimination order either way.
    for (const auto& pf : m_factors) {
        const int dimK = m_dims[(size_t)pf.k];
        for (const auto& pr : pf.L) {
            const int dimP = m_dims[(size_t)pr.first];
            y[(size_t)pr.first].head(dimP).noalias() -= pr.second.topLeftCorner(dimP, dimK) * y[(size_t)pf.k].head(dimK);
        }
    }

    // Backward substitution (U x = y), reverse elimination order: x_k = Dinv_k * (y_k - sum_q
    // U_{k,q} * x_q) -- subtract first, then apply Dinv, using the stored raw U blocks.
    std::vector<Vectorr<6>> x = y;
    for (auto it = m_factors.rbegin(); it != m_factors.rend(); ++it) {
        const int dimK = m_dims[(size_t)it->k];
        for (const auto& pr : it->U) {
            const int dimQ = m_dims[(size_t)pr.first];
            x[(size_t)it->k].head(dimK).noalias() -= pr.second.topLeftCorner(dimK, dimQ) * x[(size_t)pr.first].head(dimQ);
        }
        x[(size_t)it->k].head(dimK) = (it->Dinv.topLeftCorner(dimK, dimK) * x[(size_t)it->k].head(dimK)).eval();
    }

    VectorXr out(offset[(size_t)n]);
    for (int i = 0; i < n; ++i) out.segment(offset[(size_t)i], m_dims[(size_t)i]) = x[(size_t)i].head(m_dims[(size_t)i]);
    return out;
}

}  // namespace cardillo::misc
