#pragma once

#include <vector>
#include <string>
#include <memory>
#include <random>
#include <chrono>
#include <cstdint>
#include <iostream>

// Use the same scalar type alias as the engine when possible
#include <Eigen/Core>

namespace bench {

using Real = double;

// Simple timer utility
struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start;
    void tic() { start = Clock::now(); }
    double tocMs() const {
        auto d = std::chrono::duration<double, std::milli>(Clock::now() - start);
        return d.count();
    }
};

// Contract for block matrix/vector ops resembling contact matrix structure
// Supports mixed per-body DOF (3 for point mass, 6 for rigid). Each contact row is 3x(doi+doj).
struct IBlockOps {
    virtual ~IBlockOps() = default;
    virtual std::string name() const = 0;

    // Build a matrix with m contact rows and n bodies, each contact couples two bodies (bi, bj)
    // dofPerBody: vector of 3 or 6 per body; contacts define pairs (i,j)
    virtual void buildSystem(const std::vector<int>& dofPerBody,
                             const std::vector<std::pair<int,int>>& contacts,
                             std::uint32_t seed) = 0;

    // y = A * x (block-matrix times block-vector) where x is per-body 3-vector (e.g., velocities)
    virtual void mul(const std::vector<Real>& x, std::vector<Real>& y) const = 0;

    // z_k = row_k dot x (treating row_k as concatenation of 3x3 blocks)
    virtual Real rowDot(std::size_t k, const std::vector<Real>& x) const = 0;

    // y += A^T * w (contact-space to body-space, used for assembling forces)
    virtual void mulTransposeAcc(const std::vector<Real>& w, std::vector<Real>& y) const = 0;

    // Optionally compute A^T A x for CG-like test
    virtual void normalMatrixMul(const std::vector<Real>& x, std::vector<Real>& y) const = 0;
};

// Utility to generate random contacts over numBodies (pairs i<j)
inline std::vector<std::pair<int,int>> generateContacts(std::size_t numBodies, std::size_t numContacts, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, static_cast<int>(numBodies) - 1);
    std::vector<std::pair<int,int>> c;
    c.reserve(numContacts);
    for (std::size_t k = 0; k < numContacts; ++k) {
        int i = dist(rng), j = dist(rng);
        if (i == j) j = (j + 1) % static_cast<int>(numBodies);
        if (i > j) std::swap(i, j);
        c.emplace_back(i, j);
    }
    return c;
}

// Init vector with random values
inline void fillRandom(std::vector<Real>& v, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> U(-1.0, 1.0);
    for (auto& e : v) e = U(rng);
}

// Generate per-body DOF vector with a fixed fraction of rigids (6 DOF) vs points (3 DOF)
inline std::vector<int> generateDofPerBody(std::size_t numBodies, double rigidFrac, std::uint32_t seed) {
    std::vector<int> d(numBodies, 3);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (std::size_t i = 0; i < numBodies; ++i) {
        d[i] = (U(rng) < rigidFrac) ? 6 : 3;
    }
    return d;
}

// Compute DOF offsets and total DOF
inline std::vector<int> computeOffsets(const std::vector<int>& dofPerBody, int& totalDofOut) {
    std::vector<int> offs(dofPerBody.size(), 0);
    int acc = 0;
    for (std::size_t i = 0; i < dofPerBody.size(); ++i) {
        offs[i] = acc;
        acc += dofPerBody[i];
    }
    totalDofOut = acc;
    return offs;
}

} // namespace bench
