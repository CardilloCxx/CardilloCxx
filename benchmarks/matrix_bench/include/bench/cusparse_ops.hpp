#pragma once

#ifdef BENCH_HAVE_CUSPARSE

#include "common.hpp"
#include <cuda_runtime_api.h>
#include <cusparse.h>
#include <stdexcept>

namespace bench {

class CuSPARSEOps : public IBlockOps {
public:
    CuSPARSEOps() {
        if (cusparseCreate(&handle) != CUSPARSE_STATUS_SUCCESS) {
            throw std::runtime_error("cusparseCreate failed");
        }
    }
    ~CuSPARSEOps() override {
        release();
        if (handle) cusparseDestroy(handle);
    }

    std::string name() const override { return "cuSPARSE"; }

    void buildSystem(const std::vector<int>& dofPerBody,
                     const std::vector<std::pair<int,int>>& contacts,
                     std::uint32_t seed) override {
        (void)seed;
        dof = dofPerBody;
        rowsPairs = contacts;
        nc = (int)rowsPairs.size();
        offsets = computeOffsets(dof, ndof);

        const int m = 3 * nc;
        const int n = ndof;

        // Build CSR on host: each row has up to 2 non-zeros (i and j contributions)
        hRowPtr.assign(m + 1, 0);
        std::vector<int> cols;
        std::vector<Real> vals;
        cols.reserve((size_t)m * 2);
        vals.reserve((size_t)m * 2);
        int nnz = 0;
        for (int k = 0; k < nc; ++k) {
            int i = rowsPairs[k].first, j = rowsPairs[k].second;
            int di = dof[i], dj = dof[j];
            Real si = Real(1 + (i%5)*0.1);
            Real sj = Real(1 + (j%5)*0.1);
            for (int r = 0; r < 3; ++r) {
                // start of row
                hRowPtr[3*k + r] = nnz;
                // body i contribution
                if (r < di && r < 3) { cols.push_back(offsets[i] + r); vals.push_back(si); ++nnz; }
                // body j contribution
                if (r < dj && r < 3) { cols.push_back(offsets[j] + r); vals.push_back(-sj); ++nnz; }
            }
        }
        hRowPtr[m] = nnz;
        hColInd = std::move(cols);
        hVals = std::move(vals);

        // Allocate device memory and upload
        release();
        cudaMalloc((void**)&dRowPtr, (size_t)(m + 1) * sizeof(int));
        cudaMalloc((void**)&dColInd, (size_t)nnz * sizeof(int));
        cudaMalloc((void**)&dVals, (size_t)nnz * sizeof(Real));
        cudaMalloc((void**)&dx, (size_t)n * sizeof(Real));
        cudaMalloc((void**)&dy, (size_t)m * sizeof(Real));
        cudaMemcpy(dRowPtr, hRowPtr.data(), (size_t)(m + 1) * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dColInd, hColInd.data(), (size_t)nnz * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dVals, hVals.data(), (size_t)nnz * sizeof(Real), cudaMemcpyHostToDevice);

        // Create matrix and vector descriptors
        mRows = m; nCols = n; nnzCount = nnz;
        cusparseCreateCsr(&spMat,
                          mRows, nCols, nnzCount,
                          dRowPtr, dColInd, dVals,
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                          CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
        cusparseCreateDnVec(&vecX, nCols, dx, CUDA_R_64F);
        cusparseCreateDnVec(&vecY, mRows, dy, CUDA_R_64F);

        // Query buffer size for SpMV (NoTrans) and (Trans)
        const double alpha = 1.0, beta = 0.0;
        size_t buffer1 = 0, buffer2 = 0;
    cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                &alpha, spMat, vecX, &beta, vecY, CUDA_R_64F,
                CUSPARSE_SPMV_ALG_DEFAULT, &buffer1);
    cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_TRANSPOSE,
                &alpha, spMat, vecY, &beta, vecX, CUDA_R_64F,
                CUSPARSE_SPMV_ALG_DEFAULT, &buffer2);
        dBufferSize = std::max(buffer1, buffer2);
        cudaMalloc(&dBuffer, dBufferSize);
    }

    void mul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        // y = A*x
        y.resize(mRows);
        cudaMemcpy(dx, x.data(), (size_t)nCols * sizeof(Real), cudaMemcpyHostToDevice);
        const double alpha = 1.0, beta = 0.0;
    cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
             &alpha, spMat, vecX, &beta, vecY, CUDA_R_64F,
             CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);
        cudaMemcpy(y.data(), dy, (size_t)mRows * sizeof(Real), cudaMemcpyDeviceToHost);
    }

    Real rowDot(std::size_t k, const std::vector<Real>& x) const override {
        // Compute y = A*x once then sum the 3 entries for row k
        if (!(ycache_valid && lastXHostPtr == x.data() && lastXSize == x.size())) {
            mul(x, ycache);
            ycache_valid = true; lastXHostPtr = x.data(); lastXSize = x.size();
        }
        const size_t base = 3*k;
        return ycache[base] + ycache[base+1] + ycache[base+2];
    }

    void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const override {
        // y += A^T * w
        // Reuse descriptors by swapping vecX/vecY roles temporarily
        // Upload w into dy, compute tmp = A^T * w into dx, download and accumulate
        cudaMemcpy(dy, w.data(), (size_t)mRows * sizeof(Real), cudaMemcpyHostToDevice);
        const double alpha = 1.0, beta = 0.0;
        // Create temporary descriptors pointing to dy (input) and dx (output)
        cusparseDnVecDescr_t inW = nullptr, outG = nullptr;
        cusparseCreateDnVec(&inW, mRows, dy, CUDA_R_64F);
        cusparseCreateDnVec(&outG, nCols, dx, CUDA_R_64F);
    cusparseSpMV(handle, CUSPARSE_OPERATION_TRANSPOSE,
             &alpha, spMat, inW, &beta, outG, CUDA_R_64F,
             CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);
        std::vector<Real> ghost(nCols);
        cudaMemcpy(ghost.data(), dx, (size_t)nCols * sizeof(Real), cudaMemcpyDeviceToHost);
        if ((int)y.size() != nCols) y.resize(nCols, Real(0));
        for (int i = 0; i < nCols; ++i) y[i] += ghost[i];
        ycache_valid = false;
        cusparseDestroyDnVec(inW);
        cusparseDestroyDnVec(outG);
    }

    void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const override {
        std::vector<Real> tmp(mRows);
        mul(x, tmp);
        y.assign(nCols, Real(0));
        mulTransposeAcc(tmp, y);
    }

private:
    void release() {
        if (spMat) { cusparseDestroySpMat(spMat); spMat = nullptr; }
        if (vecX) { cusparseDestroyDnVec(vecX); vecX = nullptr; }
        if (vecY) { cusparseDestroyDnVec(vecY); vecY = nullptr; }
        if (dRowPtr) { cudaFree(dRowPtr); dRowPtr = nullptr; }
        if (dColInd) { cudaFree(dColInd); dColInd = nullptr; }
        if (dVals) { cudaFree(dVals); dVals = nullptr; }
        if (dx) { cudaFree(dx); dx = nullptr; }
        if (dy) { cudaFree(dy); dy = nullptr; }
        if (dBuffer) { cudaFree(dBuffer); dBuffer = nullptr; dBufferSize = 0; }
    }

    // sizes
    int nc = 0, ndof = 0;
    int mRows = 0, nCols = 0, nnzCount = 0;
    std::vector<int> dof, offsets;
    std::vector<std::pair<int,int>> rowsPairs;

    // host CSR
    std::vector<int> hRowPtr, hColInd;
    std::vector<Real> hVals;

    // device storage
    int* dRowPtr = nullptr;
    int* dColInd = nullptr;
    Real* dVals = nullptr;
    Real* dx = nullptr;
    Real* dy = nullptr;
    void* dBuffer = nullptr;
    size_t dBufferSize = 0;

    // cuSPARSE handles
    cusparseHandle_t handle = nullptr;
    cusparseSpMatDescr_t spMat = nullptr;
    cusparseDnVecDescr_t vecX = nullptr;
    cusparseDnVecDescr_t vecY = nullptr;

    // cache for rowDot
    mutable std::vector<Real> ycache;
    mutable bool ycache_valid = false;
    mutable const Real* lastXHostPtr = nullptr;
    mutable size_t lastXSize = 0;
};

} // namespace bench

#endif // BENCH_HAVE_CUSPARSE
