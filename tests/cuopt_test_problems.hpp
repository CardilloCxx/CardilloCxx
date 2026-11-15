#pragma once

// Cone programs min 1/2 x'Px + q'x s.t. Ax + s = b, s in K with a known KKT point, for the cuOpt
// tests: x, s and z are chosen first (s and z complementary in every cone block), then
// q = -(Px + A'z) and b = Ax + s.

#include <random>
#include <vector>

#include <Eigen/SparseCore>

#include <cardillo/physics/solver/conic/cuopt/cuopt_model.hpp>

namespace cuopt_test {

using cardillo::CscMatrix;
using cardillo::real_t;
using cardillo::Vector3r;
using cardillo::VectorXr;
using cardillo::solver::cuopt::ConeDims;

struct KnownProblem {
    CscMatrix P;       // upper triangle
    CscMatrix P_full;  // both triangles
    VectorXr q;
    CscMatrix A;
    VectorXr b;
    ConeDims dims;
    VectorXr x, s, z;  // KKT point
    real_t value{0};   // optimal value 1/2 x'Px + q'x
};

/**
 * @brief Random problem with n_coupled variables with positive curvature and n_decoupled variables
 * with zero columns of P and A (stored as explicit zeros) and no cost, like the percussions of
 * perfect constraints. P is diagonal or has 2x2 couplings. The second-order cone blocks cycle
 * through: s interior (z = 0), z interior (s = 0), both on the boundary with opposite tails, z on
 * the boundary with s = 0 (a sticking contact).
 */
inline KnownProblem randomProblem(unsigned seed, int n_coupled, int n_decoupled, const ConeDims& dims, bool diagonal) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real_t> uni(-1, 1), pos(0.5, 1.5);
    const int n = n_coupled + n_decoupled;
    const int m = dims.rows();

    KnownProblem pr;
    pr.dims = dims;

    // P = D + sum of w w' with w = a e_k + c e_{k+1}: positive definite on the coupled variables.
    std::vector<Eigen::Triplet<real_t>> tp;
    for (int j = 0; j < n_coupled; ++j) tp.emplace_back(j, j, pos(rng));
    if (!diagonal)
        for (int k = 0; k + 1 < n_coupled; k += 2) {
            const real_t a = uni(rng), c = uni(rng);
            tp.emplace_back(k, k, a * a);
            tp.emplace_back(k + 1, k + 1, c * c);
            tp.emplace_back(k, k + 1, a * c);
            tp.emplace_back(k + 1, k, a * c);
        }
    for (int j = n_coupled; j < n; ++j) tp.emplace_back(j, j, 0.0);
    pr.P_full.resize(n, n);
    pr.P_full.setFromTriplets(tp.begin(), tp.end());
    pr.P = pr.P_full.triangularView<Eigen::Upper>();

    std::vector<Eigen::Triplet<real_t>> ta;
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < 3; ++k) ta.emplace_back(i, static_cast<int>(rng() % static_cast<unsigned>(n_coupled)), uni(rng));
    for (int j = n_coupled; j < n && m > 0; ++j) ta.emplace_back(j % m, j, 0.0);
    pr.A.resize(m, n);
    pr.A.setFromTriplets(ta.begin(), ta.end());

    pr.x = VectorXr::Zero(n);
    for (int j = 0; j < n_coupled; ++j) pr.x[j] = uni(rng);
    pr.s = VectorXr::Zero(m);
    pr.z = VectorXr::Zero(m);
    int row = 0;
    for (; row < dims.zero; ++row) pr.z[row] = uni(rng);
    for (int i = 0; i < dims.nonneg; ++i, ++row) (i % 2 ? pr.s[row] : pr.z[row]) = pos(rng);
    for (std::size_t k = 0; k < dims.soc.size(); ++k) {
        const int d = dims.soc[k];
        auto sb = pr.s.segment(row, d);
        auto zb = pr.z.segment(row, d);
        VectorXr u(d - 1);
        for (int i = 0; i < d - 1; ++i) u[i] = uni(rng);
        const real_t beta = pos(rng);
        switch (d == 1 ? k % 2 : k % 4) {
            case 0:  // s interior, z = 0
                sb[0] = u.norm() + pos(rng);
                sb.tail(d - 1) = u;
                break;
            case 1:  // z interior, s = 0
                zb[0] = u.norm() + pos(rng);
                zb.tail(d - 1) = u;
                break;
            case 2:  // both on the boundary: s = (|u|, u), z = beta (|u|, -u)
                sb[0] = u.norm();
                sb.tail(d - 1) = u;
                zb[0] = beta * u.norm();
                zb.tail(d - 1) = -beta * u;
                break;
            default:  // z on the boundary, s = 0
                zb[0] = beta * u.norm();
                zb.tail(d - 1) = beta * u;
                break;
        }
        row += d;
    }

    pr.q = -(pr.P_full * pr.x + pr.A.transpose() * pr.z);
    pr.b = pr.A * pr.x + pr.s;
    pr.value = 0.5 * pr.x.dot(pr.P_full * pr.x) + pr.q.dot(pr.x);
    return pr;
}

/**
 * @brief Projection of p onto L^d as the cone program min 1/2 ||x||^2 - p'x s.t. -x + s = 0,
 * s in L^d, with the closed-form solution x = proj(p), s = x, z = x - p.
 */
inline KnownProblem coneProjection(const VectorXr& p) {
    const int d = static_cast<int>(p.size());
    KnownProblem pr;
    pr.dims.soc = {d};
    pr.P_full.resize(d, d);
    pr.P_full.setIdentity();
    pr.P = pr.P_full;
    pr.A = -pr.P_full;
    pr.q = -p;
    pr.b = VectorXr::Zero(d);

    const real_t t = p[0];
    const real_t r = d > 1 ? p.tail(d - 1).norm() : real_t(0);
    pr.x = VectorXr::Zero(d);
    if (r <= t) {
        pr.x = p;  // inside the cone
    } else if (r > -t) {
        // onto the boundary: ((t + r) / 2) (1, u / r)
        const real_t a = 0.5 * (t + r);
        pr.x[0] = a;
        if (d > 1) pr.x.tail(d - 1) = a / r * p.tail(d - 1);
    }  // else in the polar cone: x = 0
    pr.s = pr.x;
    pr.z = pr.x - p;
    pr.value = 0.5 * pr.x.squaredNorm() - p.dot(pr.x);
    return pr;
}

}  // namespace cuopt_test
