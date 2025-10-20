#pragma once

#ifdef BENCH_HAVE_CUBLAS

#include "common.hpp"
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <stdexcept>

namespace bench {

class CuBLASOps : public IBlockOps {
public:
    CuBLASOps() {
        if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cublasCreate failed");
    }
    ~CuBLASOps() override {
        release();
        if (handle) cublasDestroy(handle);
    }

    std::string name() const override { return "cuBLAS"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed;
        dof = dofPerBody;
        rows = contacts;
        nc = (int)rows.size();
        offsets = computeOffsets(dof, ndof);
        // Build host A column-major (3*nc) x ndof
        std::vector<Real> Ahost((size_t)3*nc * (size_t)ndof, Real(0));
        auto idx = [&](int r, int c){ return (size_t)c * (size_t)(3*nc) + (size_t)r; };
        for (int k = 0; k < nc; ++k) {
            int i = rows[k].first, j = rows[k].second;
            int di = dof[i], dj = dof[j];
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < std::min(3, di); ++c) if (r==c)
                    Ahost[idx(3*k + r, offsets[i] + c)] = si;
                for (int c = 0; c < std::min(3, dj); ++c) if (r==c)
                    Ahost[idx(3*k + r, offsets[j] + c)] = -sj;
            }
        }
        // Reserve device buffers
        release();
        sizeA = (size_t)3*nc * (size_t)ndof;
        cudaMalloc((void**)&dA, sizeA * sizeof(Real));
        cudaMalloc((void**)&dx, (size_t)ndof * sizeof(Real));
        cudaMalloc((void**)&dy, (size_t)3*nc * sizeof(Real));
        cudaMalloc((void**)&dg, (size_t)ndof * sizeof(Real));
        cudaMemcpy(dA, Ahost.data(), sizeA * sizeof(Real), cudaMemcpyHostToDevice);
        ycache_valid = false;
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        y.resize(3*nc);
        cudaMemcpy(dx, x.data(), (size_t)ndof * sizeof(Real), cudaMemcpyHostToDevice);
        const Real alpha = 1.0, beta = 0.0;
        cublasDgemv(handle, CUBLAS_OP_N, 3*nc, ndof, &alpha, dA, 3*nc, dx, 1, &beta, dy, 1);
        cudaMemcpy(y.data(), dy, (size_t)3*nc * sizeof(Real), cudaMemcpyDeviceToHost);
        // update caches so rowDot can reuse without recompute
        ycache = y; ycache_valid = true;
        xSynced = true; lastXHostPtr = x.data(); lastXSize = x.size();
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        // Compute y = A*x once (if cache invalid or x changed), then sum 3 entries
        if (!(ycache_valid && lastXHostPtr == x.data() && lastXSize == x.size())) {
            if (!xSynced || lastXHostPtr != x.data() || lastXSize != x.size()) {
                cudaMemcpy(dx, x.data(), (size_t)ndof * sizeof(Real), cudaMemcpyHostToDevice);
                xSynced = true; lastXHostPtr = x.data(); lastXSize = x.size();
            }
            const Real alpha = 1.0, beta = 0.0;
            cublasDgemv(handle, CUBLAS_OP_N, 3*nc, ndof, &alpha, dA, 3*nc, dx, 1, &beta, dy, 1);
            ycache.resize((size_t)3*nc);
            cudaMemcpy(ycache.data(), dy, (size_t)3*nc * sizeof(Real), cudaMemcpyDeviceToHost);
            ycache_valid = true;
        }
        const size_t base = 3*k;
        return ycache[base] + ycache[base+1] + ycache[base+2];
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        // y += A^T * w
        cudaMemcpy(dy, w.data(), (size_t)3*nc * sizeof(Real), cudaMemcpyHostToDevice);
        const Real alpha = 1.0, beta = 0.0;
        cublasDgemv(handle, CUBLAS_OP_T, 3*nc, ndof, &alpha, dA, 3*nc, dy, 1, &beta, dg, 1);
        std::vector<Real> ghost(ndof);
        cudaMemcpy(ghost.data(), dg, (size_t)ndof * sizeof(Real), cudaMemcpyDeviceToHost);
        if ((int)y.size() != ndof) y.resize(ndof, Real(0));
        for (int i = 0; i < ndof; ++i) y[i] += ghost[i];
        ycache_valid = false;
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        std::vector<Real> tmp(3*nc, Real(0));
        mul(x, tmp);
        y.assign(ndof, Real(0));
        mulTransposeAcc(tmp, y);
    }

private:
    void release() {
        if (dA) { cudaFree(dA); dA = nullptr; }
        if (dx) { cudaFree(dx); dx = nullptr; }
        if (dy) { cudaFree(dy); dy = nullptr; }
        if (dg) { cudaFree(dg); dg = nullptr; }
    }
    int nc = 0, ndof = 0;
    size_t sizeA = 0;
    std::vector<int> dof, offsets;
    std::vector<std::pair<int,int>> rows;
    cublasHandle_t handle = nullptr;
    // device buffers
    Real* dA = nullptr; Real* dx = nullptr; Real* dy = nullptr; Real* dg = nullptr;
    // cache for device x during rowDot sequences
    mutable bool xSynced = false;
    mutable const Real* lastXHostPtr = nullptr;
    mutable size_t lastXSize = 0;
    // cache for rowDot
    mutable std::vector<Real> ycache;
    mutable bool ycache_valid = false;
};

} // namespace bench

#endif // BENCH_HAVE_CUBLAS
