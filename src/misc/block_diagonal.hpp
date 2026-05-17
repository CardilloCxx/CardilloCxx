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
        n_ += block.rows();
    }

    BlockDiagonal calculateInverse() {
        BlockDiagonal inverse;

        for (const auto& block : blocks_) {
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

        return inverse;
    }

    VectorXr operator*(const VectorXr& vec) const {
        if (vec.size() != n_) throw std::invalid_argument("Dimension mismatch." + std::to_string(vec.size()) + " vs " + std::to_string(n_));

        VectorXr result(n_);

        int offset = 0;
        for (const auto& block : blocks_) {
            int blockSize = block.rows();
            result.segment(offset, blockSize) = block * vec.segment(offset, blockSize);
            offset += blockSize;
        }
        return result;
    }

   private:
    size_t n_;
    std::vector<MatrixXXr> blocks_;
};

}  // namespace cardillo