#pragma once
#include <Eigen/Cholesky>
#include <Eigen/LU>
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
 * @brief Returns a well-defined inverse of a general (not necessarily symmetric or positive
 * definite) small square block, falling back gracefully when singular: PartialPivLU -> FullPivLU
 * -> diagonal-only inverse. Never throws. Used for BlockSparseLDLT's non-symmetric (block-LU) mode
 * -- e.g. a rigid body's effective mass block under implicit gyroscopic forces, which is generally
 * non-symmetric and so cannot use invertSmallSpd()'s Cholesky/LDLT chain (both require symmetry).
 */
inline MatrixXXr invertSmallGeneral(const MatrixXXr& block) {
    MatrixXXr I = MatrixXXr::Identity(block.rows(), block.cols());
    Eigen::PartialPivLU<MatrixXXr> plu(block);
    // PartialPivLU has no .info() (it always "succeeds" numerically), so check conditioning
    // directly instead -- a tiny or non-finite determinant means the pivot was effectively
    // singular and the solve below would amplify rounding error unboundedly.
    const real_t det = plu.determinant();
    if (std::isfinite(det) && std::abs(det) > (real_t)1e-300) return plu.solve(I);

    Eigen::FullPivLU<MatrixXXr> flu(block);
    if (flu.isInvertible()) return flu.solve(I);

    std::cout << "Warning: Block singular under general (non-symmetric) inversion, using diagonal inverse only\n";
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
        VectorXr result;
        applyTo(vec, result);
        return result;
    }

    /**
     * @brief In-place variant of operator*: writes into @p out, which may alias @p in. Avoids the
     * heap allocation operator* incurs on every call, which matters when applied once per
     * CG/PJ iteration -- @p out is expected to be a buffer reused across iterations by the caller.
     * Every block in this codebase's usage is <=6x6 (contact/spring/damper rows), so the per-block
     * product is staged through a fixed Vectorr<6> stack scratch even when out aliases in.
     */
    void applyTo(const VectorXr& in, VectorXr& out) const {
        if ((std::size_t)in.size() != n_) throw std::invalid_argument("Dimension mismatch." + std::to_string(in.size()) + " vs " + std::to_string(n_));
        if (&out != &in) out.resize((Eigen::Index)n_);

        if (is_diag_ && (std::size_t)diag_.size() == n_) {
            if (&out == &in) out.array() *= diag_.array();
            else out.noalias() = diag_.cwiseProduct(in);
            return;
        }

        Vectorr<6> scratch;
        int offset = 0;
        for (const auto& block : blocks_) {
            const int blockSize = block.rows();
            scratch.head(blockSize).noalias() = block * in.segment(offset, blockSize);
            out.segment(offset, blockSize) = scratch.head(blockSize);
            offset += blockSize;
        }
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