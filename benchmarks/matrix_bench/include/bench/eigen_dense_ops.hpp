#pragma once

#include "common.hpp"
#include <Eigen/Dense>

namespace bench {

// Dense block matrix using Eigen::MatrixX<Real>, assembled sparsely then compressed
class EigenDenseOps : public IBlockOps {
public:
    std::string name() const override { return "EigenDense"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed; // deterministic assembly
        nb = static_cast<int>(dofPerBody.size());
        nc = static_cast<int>(contacts.size());
        dof = dofPerBody;
        offsets = computeOffsets(dof, ndof);
        A = Eigen::MatrixX<Real>::Zero(3*nc, ndof);
        rows.clear(); rows.reserve(nc);
        for (int k = 0; k < nc; ++k) {
            int i = contacts[k].first;
            int j = contacts[k].second;
            rows.emplace_back(i,j);
            int di = dof[i], dj = dof[j];
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            // Place simple diagonal in the first min(3,di) columns
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < std::min(3, di); ++c) {
                    if (r == c) A(3*k + r, offsets[i] + c) = si;
                }
                for (int c = 0; c < std::min(3, dj); ++c) {
                    if (r == c) A(3*k + r, offsets[j] + c) = -sj;
                }
            }
        }
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::MatrixX<Real>> xm(x.data(), ndof, 1);
        Eigen::MatrixX<Real> ym = A * xm;
        y.assign(ym.data(), ym.data() + ym.size());
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        Eigen::Map<const Eigen::MatrixX<Real>> xm(x.data(), ndof, 1);
        Eigen::MatrixX<Real> r = A.block(3*(int)k, 0, 3, ndof);
        Eigen::MatrixX<Real> val = r * xm;
        // Return sum of 3 components as a scalar row dot
        return val.sum();
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::MatrixX<Real>> wm(w.data(), 3*nc, 1);
        Eigen::Map<Eigen::MatrixX<Real>> ym(y.data(), ndof, 1);
        ym.noalias() += A.transpose() * wm;
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::MatrixX<Real>> xm(x.data(), ndof, 1);
        Eigen::MatrixX<Real> ym = A.transpose() * (A * xm);
        y.assign(ym.data(), ym.data() + ym.size());
    }

private:
    int nb = 0, nc = 0, ndof = 0;
    std::vector<int> dof, offsets;
    std::vector<std::pair<int,int>> rows;
    Eigen::MatrixX<Real> A; // size (3*nc) x ndof
};

} // namespace bench
