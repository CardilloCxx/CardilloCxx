#pragma once
#include <Eigen/Sparse>
#include <memory>
#include <stdexcept>
#include <vector>
#include "types.hpp"

namespace cardillo {

class TripletMatrix {
   public:
    using Triplet = Eigen::Triplet<real_t>;
    using ColSparse = Eigen::SparseMatrix<real_t, Eigen::ColMajor>;
    using RowSparse = Eigen::SparseMatrix<real_t, Eigen::RowMajor>;

    // Constructors
    TripletMatrix() = default;

    TripletMatrix(int rows, int cols) : n_rows_(rows), n_cols_(cols) {}

    TripletMatrix(int rows, int cols, std::shared_ptr<std::vector<Triplet>> data) : n_rows_(rows), n_cols_(cols) { blocks_.push_back({0, 0, 1.0, false, std::move(data)}); }

    // Static constructors
    static TripletMatrix zero(int rows, int cols) { return TripletMatrix(rows, cols); }

    static TripletMatrix zeroLike(const TripletMatrix& other) { return TripletMatrix(other.n_rows_, other.n_cols_); }

    static TripletMatrix fromDiag(const VectorXr& diag) {
        std::vector<Triplet> trips;
        trips.reserve((size_t)diag.size());
        for (int i = 0; i < diag.size(); ++i) {
            trips.emplace_back(i, i, diag[i]);
        }
        return TripletMatrix((int)diag.size(), (int)diag.size(), std::make_shared<std::vector<Triplet>>(std::move(trips)));
    }

    TripletMatrix zeroLike() const { return TripletMatrix(n_rows_, n_cols_); }

    // Horizontal concatenation
    TripletMatrix operator|(const TripletMatrix& other) const {
        TripletMatrix result = concat(other, true);
        result.invalidateCache();
        return result;
    }

    // Vertical concatenation
    TripletMatrix vConcat(const TripletMatrix& other) const {
        TripletMatrix result = concat(other, false);
        result.invalidateCache();
        return result;
    }

    // Scale operator: returns a scaled copy
    TripletMatrix operator*(real_t s) const {
        TripletMatrix result = *this;
        result.scale_ *= s;
        result.invalidateCache();
        return result;
    }

    // Transpose: returns a transposed copy
    TripletMatrix T() const {
        TripletMatrix result = *this;
        result.transposed_ = !result.transposed_;
        std::swap(result.n_rows_, result.n_cols_);
        std::swap(result.rowScales_, result.colScales_);
        result.invalidateCache();
        return result;
    }

    // Scale rows: returns a copy with rows scaled by the given vector
    TripletMatrix scaleRows(const VectorXr& rowScales) const {
        if (rowScales.size() != n_rows_) throw std::runtime_error("Row scale vector size mismatch");

        TripletMatrix result = *this;
        if (result.rowScales_.empty()) {
            result.rowScales_.assign((size_t)rowScales.size(), (real_t)1);
        }
        for (int i = 0; i < rowScales.size(); ++i) {
            result.rowScales_[(size_t)i] *= rowScales[i];
        }
        result.invalidateCache();

        return result;
    }

    // Convert to Eigen SparseMatrix with internal caching.
    const ColSparse& asSparse() const {
        if (!sparse_cache_col_) {
            auto outTriplets = materializeTriplets_();
            ColSparse mat(n_rows_, n_cols_);
            mat.setFromTriplets(outTriplets.begin(), outTriplets.end());
            mat.makeCompressed();
            sparse_cache_col_ = std::make_shared<ColSparse>(std::move(mat));
        }
        return *sparse_cache_col_;
    }

    const RowSparse& asSparseRowMajor() const {
        if (!sparse_cache_row_) {
            auto outTriplets = materializeTriplets_();
            RowSparse mat(n_rows_, n_cols_);
            mat.setFromTriplets(outTriplets.begin(), outTriplets.end());
            mat.makeCompressed();
            sparse_cache_row_ = std::make_shared<RowSparse>(std::move(mat));
        }
        return *sparse_cache_row_;
    }

    // Accessors
    int nRows() const { return n_rows_; }
    int nCols() const { return n_cols_; }

   private:
    struct Block {
        int row_offset = 0;
        int col_offset = 0;
        real_t scale = 1.0;
        bool transposed = false;
        std::shared_ptr<std::vector<Triplet>> data;
    };

    int n_rows_ = 0;
    int n_cols_ = 0;
    real_t scale_ = 1.0;
    bool transposed_ = false;
    std::vector<real_t> rowScales_;
    std::vector<real_t> colScales_;

    std::vector<Block> blocks_;

    mutable std::shared_ptr<ColSparse> sparse_cache_col_;
    mutable std::shared_ptr<RowSparse> sparse_cache_row_;

    std::vector<Triplet> materializeTriplets_() const {
        std::vector<Triplet> outTriplets;
        size_t totalSize = 0;
        for (const auto& b : blocks_) {
            totalSize += b.data->size();
        }
        outTriplets.reserve(totalSize);

        for (const auto& b : blocks_) {
            bool effective_transpose = transposed_ ^ b.transposed;
            for (const auto& t : *(b.data)) {
                int r = effective_transpose ? t.col() : t.row();
                int c = effective_transpose ? t.row() : t.col();
                const int rr = r + b.row_offset;
                const int cc = c + b.col_offset;

                real_t v = t.value() * b.scale * scale_;
                if (!rowScales_.empty()) v *= rowScales_[(size_t)rr];
                if (!colScales_.empty()) v *= colScales_[(size_t)cc];

                outTriplets.emplace_back(rr, cc, v);
            }
        }
        return outTriplets;
    }

    // Shared concat logic
    TripletMatrix concat(const TripletMatrix& other, bool horizontal) const {
        if (horizontal && n_rows_ != other.n_rows_) throw std::runtime_error("Row mismatch for horizontal concat");
        if (!horizontal && n_cols_ != other.n_cols_) throw std::runtime_error("Column mismatch for vertical concat");

        TripletMatrix result;
        result.n_rows_ = horizontal ? n_rows_ : n_rows_ + other.n_rows_;
        result.n_cols_ = horizontal ? n_cols_ + other.n_cols_ : n_cols_;

        // Fold pending matrix-level scale/transpose into block metadata so
        // concatenation preserves transforms from both operands.
        result.blocks_ = transformedBlocks_();
        auto otherBlocks = other.transformedBlocks_();
        for (const auto& b : otherBlocks) {
            Block shifted = b;
            if (horizontal)
                shifted.col_offset += n_cols_;
            else
                shifted.row_offset += n_rows_;
            result.blocks_.push_back(shifted);
        }

        result.scale_ = 1.0;
        result.transposed_ = false;

        return result;
    }

    std::vector<Block> transformedBlocks_() const {
        std::vector<Block> transformed;
        transformed.reserve(blocks_.size());
        for (const auto& b : blocks_) {
            Block out = b;
            out.scale *= scale_;
            if (transposed_) {
                out.transposed = !out.transposed;
                std::swap(out.row_offset, out.col_offset);
            }

            // Fold pending row/column scaling into block values so concat preserves semantics.
            if (!rowScales_.empty() || !colScales_.empty()) {
                auto scaledData = std::make_shared<std::vector<Triplet>>();
                scaledData->reserve(out.data->size());

                for (const auto& t : *out.data) {
                    const int lr = out.transposed ? t.col() : t.row();
                    const int lc = out.transposed ? t.row() : t.col();
                    const int rr = lr + out.row_offset;
                    const int cc = lc + out.col_offset;

                    real_t v = t.value();
                    if (!rowScales_.empty()) v *= rowScales_[(size_t)rr];
                    if (!colScales_.empty()) v *= colScales_[(size_t)cc];
                    scaledData->emplace_back(t.row(), t.col(), v);
                }

                out.data = std::move(scaledData);
            }
            transformed.push_back(out);
        }
        return transformed;
    }

    // Invalidate cached sparse matrix
    void invalidateCache() const {
        sparse_cache_col_.reset();
        sparse_cache_row_.reset();
    }
};

}  // namespace cardillo