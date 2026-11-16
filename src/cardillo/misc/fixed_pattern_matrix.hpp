#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include <Eigen/SparseCore>

#include "types.hpp"

namespace cardillo::misc {

// Sparse matrix assembled from the same sequence of entries every step, e.g. the spring and damper
// rows of DynamicsAssembler or the saddle-point matrix S of projected Jacobi. Each step adds its
// entries with add() between begin() and finish(). The sparsity pattern holds every added entry
// that has been nonzero; entries that became zero again stay in it. finish() writes the values into
// the compressed matrix through the value positions found when the pattern was built, so the
// pattern only changes when an entry becomes nonzero for the first time (it grows) or when entries
// are no longer added (bodies, constraints or used rows removed; they leave it). Entries that are
// always zero never enter it. Entries added twice at the same position are summed.
class FixedPatternMatrix {
   public:
    // Starts the entries of a step; `cols` is the number of columns of the matrix.
    void begin(Eigen::Index cols) {
        m_cols = cols;
        m_stepKeys.clear();
        m_stepValues.clear();
    }

    void add(Eigen::Index row, Eigen::Index col, real_t value) {
        m_stepKeys.push_back((std::int64_t)row * m_cols + col);
        m_stepValues.push_back(value);
    }

    // Writes the values of the step into the matrix (rows x cols); returns true if the pattern changed.
    bool finish(Eigen::Index rows) {
        const bool sameEntries = m_matrix.rows() == rows && m_matrix.cols() == m_cols && m_stepKeys == m_keys;
        if (sameEntries && writeValues()) return false;
        buildPattern(rows);
        writeValues();
        return true;
    }

    const CscMatrix& matrix() const { return m_matrix; }

   private:
    // false if a nonzero value lies outside the pattern
    bool writeValues() {
        real_t* values = m_matrix.valuePtr();
        std::fill(values, values + m_matrix.nonZeros(), (real_t)0);
        for (std::size_t e = 0; e < m_stepValues.size(); ++e) {
            const Eigen::Index position = m_positions[e];
            if (position >= 0)
                values[position] += m_stepValues[e];
            else if (m_stepValues[e] != (real_t)0)
                return false;
        }
        return true;
    }

    // The added entries that are nonzero now or in the current pattern.
    void buildPattern(Eigen::Index rows) {
        const CscMatrix previous = std::move(m_matrix);
        auto inPrevious = [&previous](CscMatrix::StorageIndex row, Eigen::Index col) {
            if (col >= previous.cols() || row >= previous.rows()) return false;
            const auto* inner = previous.innerIndexPtr();
            const auto* end = inner + previous.outerIndexPtr()[col + 1];
            const auto* found = std::lower_bound(inner + previous.outerIndexPtr()[col], end, row);
            return found != end && *found == row;
        };

        std::vector<Eigen::Triplet<real_t>> triplets;
        triplets.reserve(m_stepKeys.size());
        for (std::size_t e = 0; e < m_stepKeys.size(); ++e) {
            const auto row = (CscMatrix::StorageIndex)(m_stepKeys[e] / m_cols);
            const auto col = (Eigen::Index)(m_stepKeys[e] % m_cols);
            if (m_stepValues[e] != (real_t)0 || inPrevious(row, col)) triplets.emplace_back((int)row, (int)col, (real_t)0);
        }
        m_matrix = CscMatrix(rows, m_cols);
        m_matrix.setFromTriplets(triplets.begin(), triplets.end());
        m_matrix.makeCompressed();

        // Value position of every entry: its row in the sorted row indices of its column, -1 if absent.
        const auto* outer = m_matrix.outerIndexPtr();
        const auto* inner = m_matrix.innerIndexPtr();
        m_positions.resize(m_stepKeys.size());
        for (std::size_t e = 0; e < m_stepKeys.size(); ++e) {
            const auto row = (CscMatrix::StorageIndex)(m_stepKeys[e] / m_cols);
            const auto col = (Eigen::Index)(m_stepKeys[e] % m_cols);
            const auto* end = inner + outer[col + 1];
            const auto* found = std::lower_bound(inner + outer[col], end, row);
            m_positions[e] = (found != end && *found == row) ? found - inner : -1;
        }
        m_keys = m_stepKeys;
    }

    Eigen::Index m_cols{0};
    CscMatrix m_matrix;
    std::vector<std::int64_t> m_keys;       // sequence of entries (row * cols + col) of the pattern
    std::vector<Eigen::Index> m_positions;  // entry -> index into m_matrix.valuePtr(), -1: not in the pattern
    std::vector<std::int64_t> m_stepKeys;   // entries of the current step
    std::vector<real_t> m_stepValues;
};

}  // namespace cardillo::misc
