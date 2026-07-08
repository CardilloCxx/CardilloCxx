#include "block_sparse_ldlt.hpp"

#include <limits>
#include <set>

#include "block_diagonal.hpp"  // invertSmallSpd

namespace cardillo::misc {

void BlockSparseLDLT::build(std::vector<int> dims, std::vector<MatrixXXr> diagBlocks, const std::vector<std::array<int, 2>>& edgeNodes, const std::vector<MatrixXXr>& edgeBlocks) {
    m_dims = std::move(dims);
    m_diag = std::move(diagBlocks);
    m_off.assign(m_dims.size(), {});
    m_factors.clear();

    for (size_t e = 0; e < edgeNodes.size(); ++e) {
        const int i = edgeNodes[e][0];
        const int j = edgeNodes[e][1];
        m_off[(size_t)i][j] = edgeBlocks[e];
        m_off[(size_t)j][i] = edgeBlocks[e].transpose();
    }
}

int BlockSparseLDLT::totalDim() const {
    int total = 0;
    for (int d : m_dims) total += d;
    return total;
}

std::vector<int> BlockSparseLDLT::computeGreedyMinDegreeOrder() const {
    const int n = (int)m_dims.size();
    std::vector<std::set<int>> adj(n);
    for (int i = 0; i < n; ++i)
        for (const auto& kv : m_off[(size_t)i]) adj[(size_t)i].insert(kv.first);

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
    m_factors.clear();
    m_factors.reserve(m_dims.size());
    std::vector<bool> eliminated(m_dims.size(), false);

    for (int k : m_order) {
        PivotFactor pf;
        pf.k = k;
        pf.Dinv = invertSmallSpd(m_diag[(size_t)k]);

        // Still-active neighbors of k (copied out before mutating anything -- the update loop
        // below writes into m_off[p]/m_off[q] for p,q in this list, never into m_off[k] itself).
        std::vector<int> nbrs;
        nbrs.reserve(m_off[(size_t)k].size());
        for (const auto& kv : m_off[(size_t)k]) {
            if (!eliminated[(size_t)kv.first]) nbrs.push_back(kv.first);
        }

        pf.L.reserve(nbrs.size());
        for (int p : nbrs) {
            MatrixXXr Lpk = m_off[(size_t)p].at(k) * pf.Dinv;  // dims[p] x dims[k]
            pf.L.emplace_back(p, std::move(Lpk));
        }

        // Right-looking Schur update over every pair (p, q) of still-active neighbors, including
        // p==q (the diagonal update) -- this is the block generalization of one elimination step
        // of a scalar sparse Cholesky/Thomas algorithm. Creates fill-in edges for p!=q pairs not
        // already connected.
        for (size_t ii = 0; ii < nbrs.size(); ++ii) {
            const int p = nbrs[ii];
            const MatrixXXr& Lp = pf.L[ii].second;  // dims[p] x dims[k]
            for (size_t jj = ii; jj < nbrs.size(); ++jj) {
                const int q = nbrs[jj];
                const MatrixXXr& Bkq = m_off[(size_t)k].at(q);  // dims[k] x dims[q]
                MatrixXXr update = Lp * Bkq;                    // dims[p] x dims[q]

                if (p == q) {
                    m_diag[(size_t)p] -= update;
                    m_diag[(size_t)p] = (real_t)0.5 * (m_diag[(size_t)p] + m_diag[(size_t)p].transpose());  // keep numerically symmetric
                } else {
                    auto it = m_off[(size_t)p].find(q);
                    if (it != m_off[(size_t)p].end()) {
                        it->second -= update;
                    } else {
                        m_off[(size_t)p][q] = -update;  // fill-in
                    }
                    m_off[(size_t)q][p] = m_off[(size_t)p].at(q).transpose();
                }
            }
        }

        eliminated[(size_t)k] = true;
        m_factors.push_back(std::move(pf));
    }
}

VectorXr BlockSparseLDLT::solve(const VectorXr& rhs) const {
    const int n = (int)m_dims.size();
    std::vector<int> offset(n + 1, 0);
    for (int i = 0; i < n; ++i) offset[(size_t)i + 1] = offset[(size_t)i] + m_dims[(size_t)i];

    std::vector<VectorXr> y(n);
    for (int i = 0; i < n; ++i) y[(size_t)i] = rhs.segment(offset[(size_t)i], m_dims[(size_t)i]);

    // Forward substitution (L y = rhs, L unit lower-triangular in elimination order): visiting
    // pivots in elimination order guarantees y[pf.k] already holds every correction from earlier
    // pivots before it propagates its own effect onto its (later-eliminated) neighbors.
    for (const auto& pf : m_factors) {
        for (const auto& pr : pf.L) {
            y[(size_t)pr.first].noalias() -= pr.second * y[(size_t)pf.k];
        }
    }

    std::vector<VectorXr> x = y;
    for (const auto& pf : m_factors) x[(size_t)pf.k] = pf.Dinv * y[(size_t)pf.k];

    // Backward substitution (L^T x = D^-1 y), reverse elimination order.
    for (auto it = m_factors.rbegin(); it != m_factors.rend(); ++it) {
        for (const auto& pr : it->L) {
            x[(size_t)it->k].noalias() -= pr.second.transpose() * x[(size_t)pr.first];
        }
    }

    VectorXr out(offset[(size_t)n]);
    for (int i = 0; i < n; ++i) out.segment(offset[(size_t)i], m_dims[(size_t)i]) = x[(size_t)i];
    return out;
}

}  // namespace cardillo::misc
