#pragma once
#include <Eigen/Cholesky>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "types.hpp"

namespace cardillo {

/**
 * @brief Returns a well-defined inverse of a nominally SPD square block, falling back gracefully
 * when it isn't: Cholesky (LLT) -> LDLT -> diagonal-only inverse. Never throws.
 */
inline MatrixXXr invertSmallSpd(const MatrixXXr& block) {
    MatrixXXr I = MatrixXXr::Identity(block.rows(), block.cols());
    Eigen::LLT<MatrixXXr> llt(block);
    if (llt.info() == Eigen::Success) return llt.solve(I);

    Eigen::LDLT<MatrixXXr> ldlt(block);
    if (ldlt.info() == Eigen::Success) return ldlt.solve(I);

    std::cout << "Warning: Block not positive definite, using diagonal inverse only\n";
    MatrixXXr diagInv = MatrixXXr::Zero(block.rows(), block.cols());
    diagInv.diagonal() = block.diagonal().cwiseInverse();
    return diagInv;
}

/**
 * @brief A block-diagonal matrix, stored as a list of dense square blocks.
 *
 * Used to represent the (block-diagonal) Delassus-type operators assembled by the PGS/CG
 * solvers, where each block corresponds to one constraint's local compliance/mass contribution.
 * If every block added so far is itself diagonal, `is_diag_` stays true and operations fall back
 * to a cheap element-wise diagonal path.
 */
class BlockDiagonal {
   public:
    BlockDiagonal() = default;

    /// Appends a square block; must be called in the same order as the corresponding DOF layout.
    void addBlock(const MatrixXXr& block) {
        if (block.rows() != block.cols()) throw std::invalid_argument("Block must be square.");
        blocks_.push_back(block);
        is_diag_ = is_diag_ && block.isDiagonal();
        n_ += block.rows();
    }

    /// Appends a diagonal block built from `diag`.
    void addBlockDiag(const VectorXr& diag) {
        MatrixXXr block = diag.asDiagonal();
        addBlock(block);
    }

    /**
     * @brief Returns the block-wise inverse of this matrix.
     *
     * Each block is inverted independently via Cholesky (LLT), falling back to LDLT and then to
     * a diagonal-only inverse if the block is not (numerically) positive definite.
     */
    BlockDiagonal calculateInverse() const {
        BlockDiagonal inverse;

        for (const auto& block : blocks_) {
            if (is_diag_) {
                inverse.addBlockDiag(block.diagonal().cwiseInverse());
                continue;
            }
            inverse.addBlock(invertSmallSpd(block));
        }
        inverse.calcDiag();
        return inverse;
    }

    VectorXr operator*(const VectorXr& vec) const {
        if (vec.size() != n_) throw std::invalid_argument("Dimension mismatch." + std::to_string(vec.size()) + " vs " + std::to_string(n_));

        if (is_diag_ && vec.size() == diag_.size()) {
            return diag_.cwiseProduct(vec);
        }

        VectorXr result(n_);

        int offset = 0;
        for (const auto& block : blocks_) {
            int blockSize = block.rows();
            result.segment(offset, blockSize) = block * vec.segment(offset, blockSize);
            offset += blockSize;
        }
        return result;
    }

    const std::vector<MatrixXXr>& blocks() const { return blocks_; }
    void calcDiag() {
        diag_ = VectorXr::Zero(n_);
        int offset = 0;
        for (const auto& block : blocks_) {
            diag_.segment(offset, block.rows()) = block.diagonal();
            offset += block.rows();
        }
    }
    size_t size() const { return n_; }

   private:
    size_t n_ = 0;
    std::vector<MatrixXXr> blocks_;
    bool is_diag_ = true;
    VectorXr diag_;
};

}  // namespace cardillo