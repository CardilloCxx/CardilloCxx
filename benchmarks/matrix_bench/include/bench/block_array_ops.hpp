#pragma once

#include "common.hpp"
#include <array>

namespace bench {

// Block-array of 3x(doi) per contact row: stores bi and bj indices and two dense 3x(doi) blocks (diagonal-ish)
class BlockArrayOps : public IBlockOps {
public:
    std::string name() const override { return "BlockArray"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed;
        nb = static_cast<int>(dofPerBody.size());
        rows = contacts;
        nc = static_cast<int>(rows.size());
        dof = dofPerBody;
        offsets = computeOffsets(dof, ndof);
        Bi.clear(); Bj.clear(); Bi.resize(nc); Bj.resize(nc);
        for (int k = 0; k < nc; ++k) {
            int i = rows[k].first;
            int j = rows[k].second;
            int di = dof[i], dj = dof[j];
            Bi[k].assign(3*di, Real(0));
            Bj[k].assign(3*dj, Real(0));
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            // Fill a simple diagonal across the first min(3,di) columns to mimic coupling; same for j
            int wi = std::min(3, di);
            int wj = std::min(3, dj);
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < wi; ++c) {
                    if (r == c) Bi[k][r*di + c] = si;
                }
                for (int c = 0; c < wj; ++c) {
                    if (r == c) Bj[k][r*dj + c] = -sj;
                }
            }
        }
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        y.assign(3*nc, Real(0));
        for (int k = 0; k < nc; ++k) {
            int i = rows[k].first;
            int j = rows[k].second;
            const Real* Xi = &x[offsets[i]];
            const Real* Xj = &x[offsets[j]];
            Real* Yk = &y[3*k];
            // y_k = Bi * x_i + Bj * x_j; Bi is 3 x di, Bj is 3 x dj
            int di = dof[i], dj = dof[j];
            for (int r = 0; r < 3; ++r) {
                Real acc = 0;
                for (int c = 0; c < di; ++c) acc += Bi[k][r*di + c] * Xi[c];
                for (int c = 0; c < dj; ++c) acc += Bj[k][r*dj + c] * Xj[c];
                Yk[r] = acc;
            }
        }
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        int i = rows[k].first;
        int j = rows[k].second;
        const Real* Xi = &x[offsets[i]];
        const Real* Xj = &x[offsets[j]];
        int di = dof[i], dj = dof[j];
        Real sum3 = 0;
        for (int r = 0; r < 3; ++r) {
            Real acc = 0;
            for (int c = 0; c < di; ++c) acc += Bi[k][r*di + c] * Xi[c];
            for (int c = 0; c < dj; ++c) acc += Bj[k][r*dj + c] * Xj[c];
            sum3 += acc;
        }
        return sum3;
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        // y += A^T * w
        for (int k = 0; k < nc; ++k) {
            int i = rows[k].first;
            int j = rows[k].second;
            const Real* Wk = &w[3*k];
            Real* Yi = &y[offsets[i]];
            Real* Yj = &y[offsets[j]];
            int di = dof[i], dj = dof[j];
            // Yi += Bi^T * Wk; Yj += Bj^T * Wk
            for (int c = 0; c < di; ++c) {
                Real acc_i = 0;
                for (int r = 0; r < 3; ++r) acc_i += Bi[k][r*di + c] * Wk[r];
                Yi[c] += acc_i;
            }
            for (int c = 0; c < dj; ++c) {
                Real acc_j = 0;
                for (int r = 0; r < 3; ++r) acc_j += Bj[k][r*dj + c] * Wk[r];
                Yj[c] += acc_j;
            }
        }
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        // y = A^T A x
        std::vector<Real> tmp(3*nc, Real(0));
        mul(x, tmp);
        y.assign(ndof, Real(0));
        mulTransposeAcc(tmp, y);
    }

private:
    int nb = 0, nc = 0;
    int ndof = 0;
    std::vector<std::pair<int,int>> rows;
    std::vector<int> dof;
    std::vector<int> offsets;
    // Store each 3xdi block row-major in a flat vector sized 3*di (same for dj)
    std::vector<std::vector<Real>> Bi, Bj; // per-contact 3x(doi) and 3x(doj) blocks
};

} // namespace bench
