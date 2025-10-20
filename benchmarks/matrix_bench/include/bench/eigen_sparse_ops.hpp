#pragma once

#include "common.hpp"
#include <Eigen/Sparse>

namespace bench {

class EigenSparseOps : public IBlockOps {
public:
    std::string name() const override { return "EigenSparse"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed;
        nb = static_cast<int>(dofPerBody.size());
        nc = static_cast<int>(contacts.size());
        dof = dofPerBody;
        offsets = computeOffsets(dof, ndof);
        std::vector<Eigen::Triplet<Real>> triplets;
        triplets.reserve(nc * 2 * 3 * 3); // upper bound
        rows.clear(); rows.reserve(nc);
        for (int k = 0; k < nc; ++k) {
            int i = contacts[k].first;
            int j = contacts[k].second;
            rows.emplace_back(i,j);
            int di = dof[i], dj = dof[j];
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < std::min(3, di); ++c) if (r==c)
                    triplets.emplace_back(3*k + r, offsets[i] + c, si);
                for (int c = 0; c < std::min(3, dj); ++c) if (r==c)
                    triplets.emplace_back(3*k + r, offsets[j] + c, -sj);
            }
        }
        A.resize(3*nc, ndof);
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> xm(x.data(), ndof);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> ym(3*nc);
        ym = A * xm;
        y.assign(ym.data(), ym.data() + ym.size());
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> xm(x.data(), ndof);
        // Efficient row * vector multiply via innerIterator on row k
        Real acc = 0;
        for (typename Eigen::SparseMatrix<Real>::InnerIterator it(A, (int)k); it; ++it) {
            acc += it.value() * xm(it.col());
        }
        return acc;
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> wm(w.data(), 3*nc);
        Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>> ym(y.data(), ndof);
        ym.noalias() += A.transpose() * wm;
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, 1>> xm(x.data(), ndof);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> ym(ndof);
        ym = A.transpose() * (A * xm);
        y.assign(ym.data(), ym.data() + ym.size());
    }

private:
    int nb = 0, nc = 0, ndof = 0;
    std::vector<int> dof, offsets;
    std::vector<std::pair<int,int>> rows;
    Eigen::SparseMatrix<Real> A; // (3*nc) x ndof
};

} // namespace bench
