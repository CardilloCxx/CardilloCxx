#pragma once
#include <iostream>
#include <stdexcept>
#include <vector>
#include "types.hpp"

namespace cardillo {

class BlockDiagonal {
   public:
    BlockDiagonal() = default;

    void addBlock(const MatrixXXr& block) {
        if (block.rows() != block.cols()) throw std::invalid_argument("Block must be square.");
        blocks_.push_back(block);
        is_diag_ = is_diag_ && block.isDiagonal();
        n_ += block.rows();
    }

    void addBlockDiag(const VectorXr& diag) {
        MatrixXXr block = diag.asDiagonal();
        addBlock(block);
    }

    BlockDiagonal calculateInverse() {
        BlockDiagonal inverse;

        for (const auto& block : blocks_) {
            if (is_diag_) {
                inverse.addBlockDiag(block.diagonal().cwiseInverse());
                continue;
            }

            MatrixXXr I = MatrixXXr::Identity(block.rows(), block.cols());
            Eigen::LLT<MatrixXXr> llt(block);
            if (llt.info() == Eigen::Success) {
                inverse.addBlock(llt.solve(I));
            } else {
                Eigen::LDLT<MatrixXXr> ldlt(block);
                if (ldlt.info() == Eigen::Success) {
                    inverse.addBlock(ldlt.solve(I));
                } else {
                    std::cout << "Warning: Block not positive definite, using diagonal inverse only\n";
                    MatrixXXr diagInv = MatrixXXr::Zero(block.rows(), block.cols());
                    diagInv.diagonal() = block.diagonal().cwiseInverse();
                    inverse.addBlock(diagInv);
                }
            }
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
    const void calcDiag() {
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