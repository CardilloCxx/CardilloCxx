#pragma once

#include "common.hpp"
#include <cblas.h>

namespace bench {

// CPU BLAS backend using column-major A and cblas_dgemv
class BlasOps : public IBlockOps {
public:
    std::string name() const override { return "BLAS"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed;
        dof = dofPerBody;
        rows = contacts;
        nc = (int)rows.size();
        offsets = computeOffsets(dof, ndof);
        // Build column-major dense A of size (3*nc) x ndof, allocate once
        A.resize((size_t)3*nc * (size_t)ndof, Real(0));
        auto idx = [&](int r, int c){ return (size_t)c * (size_t)(3*nc) + (size_t)r; }; // col-major
        for (int k = 0; k < nc; ++k) {
            int i = rows[k].first, j = rows[k].second;
            int di = dof[i], dj = dof[j];
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < std::min(3, di); ++c) if (r==c)
                    A[idx(3*k + r, offsets[i] + c)] = si;
                for (int c = 0; c < std::min(3, dj); ++c) if (r==c)
                    A[idx(3*k + r, offsets[j] + c)] = -sj;
            }
        }
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        y.resize(3*nc);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3*nc, ndof, 1.0,
                    A.data(), 3*nc, x.data(), 1, 0.0, y.data(), 1);
        // update cache for rowDot sweeps
        ycache = y; ycache_valid = true; lastXHostPtr = x.data(); lastXSize = x.size();
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        // Compute y = A*x once per x, then sum the 3 entries for row k
        if (!(ycache_valid && lastXHostPtr == x.data() && lastXSize == x.size())) {
            ycache.resize(3*nc);
            cblas_dgemv(CblasColMajor, CblasNoTrans, 3*nc, ndof, 1.0,
                        A.data(), 3*nc, x.data(), 1, 0.0, ycache.data(), 1);
            ycache_valid = true; lastXHostPtr = x.data(); lastXSize = x.size();
        }
        size_t base = 3*k;
        return ycache[base] + ycache[base+1] + ycache[base+2];
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        // y += A^T * w
        std::vector<Real> tmp(ndof, Real(0));
        cblas_dgemv(CblasColMajor, CblasTrans, 3*nc, ndof, 1.0,
                    A.data(), 3*nc, w.data(), 1, 0.0, tmp.data(), 1);
        if ((int)y.size() != ndof) y.resize(ndof, Real(0));
        for (int i = 0; i < ndof; ++i) y[i] += tmp[i];
        ycache_valid = false;
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        std::vector<Real> tmp(3*nc, Real(0));
        mul(x, tmp);
        y.assign(ndof, Real(0));
        mulTransposeAcc(tmp, y);
    }

private:
    int nc = 0, ndof = 0;
    std::vector<int> dof, offsets;
    std::vector<std::pair<int,int>> rows;
    std::vector<Real> A; // column-major (3*nc) x ndof
    // cache for rowDot sweeps
    mutable std::vector<Real> ycache;
    mutable bool ycache_valid = false;
    mutable const Real* lastXHostPtr = nullptr;
    mutable size_t lastXSize = 0;
};

} // namespace bench
