#pragma once

// Shared helpers for the conic-solver tests/benchmarks: a solver-independent problem description in
// the stacked form used by ClarabelAssembler (min 1/2 x'Px + q'x, A x + s = b, s in
// {0}^p x R+^l x Q^q1 x ...), reference solves with Clarabel and QOCO on exactly the same data,
// and KKT residuals.

#include <qoco.h>
#include <Eigen/SparseCore>
#include <algorithm>
#include <chrono>
#include <clarabel.hpp>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "physics/solver/mosek_backend.hpp"

namespace conic_test {

using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Vec = Eigen::VectorXd;
using cardillo::solver::mosek::ConeDims;

struct Problem {
    SpMat P;  // full symmetric
    Vec q;
    SpMat A;  // stacked [A_eq; G]
    Vec b;    // stacked [b; h]
    ConeDims dims;
};

struct Result {
    bool ok{false};
    Vec x, z;  // z: stacked [z_eq; z_cone], sign convention P x + q + A^T z = 0
    double obj{0};
    double seconds{0};
};

inline double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

inline std::vector<clarabel::SupportedConeT<double>> clarabelCones(const ConeDims& d) {
    std::vector<clarabel::SupportedConeT<double>> cones;
    if (d.zero > 0) cones.push_back(clarabel::ZeroConeT<double>(d.zero));
    if (d.nonneg > 0) cones.push_back(clarabel::NonnegativeConeT<double>(d.nonneg));
    for (int k : d.soc) cones.push_back(clarabel::SecondOrderConeT<double>(static_cast<uintptr_t>(k)));
    return cones;
}

inline Result solveClarabel(const Problem& pr) {
    const double t0 = now();
    auto settings = clarabel::DefaultSettings<double>::default_settings();
    settings.verbose = false;
    settings.tol_gap_abs = 1e-9;
    settings.tol_gap_rel = 1e-9;
    settings.tol_feas = 1e-9;
    SpMat Pu = pr.P.triangularView<Eigen::Upper>();
    Pu.makeCompressed();
    Vec q = pr.q, b = pr.b;
    clarabel::DefaultSolver<double> solver(Pu, q, pr.A, b, clarabelCones(pr.dims), settings);
    solver.solve();
    auto sol = solver.solution();
    Result r;
    r.ok = sol.status == clarabel::SolverStatus::Solved || sol.status == clarabel::SolverStatus::AlmostSolved;
    r.x = sol.x;
    r.z = sol.z;
    r.obj = sol.obj_val;
    r.seconds = now() - t0;
    return r;
}

// QOCO wants A/b (equalities) and G/h (cones) separately and P as upper triangle.
inline Result solveQoco(const Problem& pr) {
    const double t0 = now();
    const int n = static_cast<int>(pr.P.cols());
    const int p = pr.dims.zero;
    const int m = static_cast<int>(pr.A.rows()) - p;
    SpMat Pu = pr.P.triangularView<Eigen::Upper>();
    SpMat Aeq = pr.A.topRows(p);
    SpMat G = pr.A.bottomRows(m);
    Pu.makeCompressed();
    Aeq.makeCompressed();
    G.makeCompressed();
    Vec c = pr.q, b = pr.b.head(p), h = pr.b.tail(m);
    std::vector<QOCOInt> qv(pr.dims.soc.begin(), pr.dims.soc.end());

    QOCOCscMatrix Pq, Aq, Gq;
    qoco_set_csc(&Pq, n, n, Pu.nonZeros(), Pu.valuePtr(), Pu.outerIndexPtr(), Pu.innerIndexPtr());
    qoco_set_csc(&Aq, p, n, Aeq.nonZeros(), Aeq.valuePtr(), Aeq.outerIndexPtr(), Aeq.innerIndexPtr());
    qoco_set_csc(&Gq, m, n, G.nonZeros(), G.valuePtr(), G.outerIndexPtr(), G.innerIndexPtr());

    QOCOSettings settings;
    set_default_settings(&settings);
    settings.verbose = 0;
    settings.abstol = 1e-9;
    settings.reltol = 1e-9;
    QOCOSolver* solver = static_cast<QOCOSolver*>(malloc(sizeof(QOCOSolver)));
    Result r;
    if (qoco_setup(solver, n, m, p, Pu.nonZeros() ? &Pq : nullptr, c.data(), p ? &Aq : nullptr, p ? b.data() : nullptr, m ? &Gq : nullptr, m ? h.data() : nullptr, pr.dims.nonneg,
                   static_cast<QOCOInt>(qv.size()), qv.empty() ? nullptr : qv.data(), &settings) != 0) {
        free(solver);
        return r;
    }
    const QOCOInt exit = qoco_solve(solver);
    r.ok = exit == QOCO_SOLVED || exit == QOCO_SOLVED_INACCURATE;
    r.x = Eigen::Map<Vec>(solver->sol->x, n);
    r.z.resize(p + m);
    // QOCO: P x + c + A^T y + G^T z = 0  (same convention, y first)
    if (p) r.z.head(p) = Eigen::Map<Vec>(solver->sol->y, p);
    if (m) r.z.tail(m) = Eigen::Map<Vec>(solver->sol->z, m);
    r.obj = solver->sol->obj;
    qoco_cleanup(solver);
    r.seconds = now() - t0;
    return r;
}

inline Result solveMosek(cardillo::solver::mosek::ConicSolver& solver, const Problem& pr) {
    const double t0 = now();
    Result r;
    const auto status = solver.solve(pr.P, pr.q, pr.A, pr.b, pr.dims);
    r.ok = cardillo::solver::mosek::isAcceptable(status);
    if (!r.ok) return r;
    const int p = pr.dims.zero;
    r.x = solver.x();
    r.z.resize(pr.A.rows());
    if (solver.equalityDuals().size() == p)
        r.z.head(p) = solver.equalityDuals();
    else
        r.z.head(p).setConstant(std::numeric_limits<double>::quiet_NaN());  // Options::equality_duals off
    r.z.tail(pr.A.rows() - p) = solver.coneDuals();
    r.obj = solver.objective();
    r.seconds = now() - t0;
    return r;
}

inline Result solveMosek(const Problem& pr) {
    cardillo::solver::mosek::ConicSolver solver;
    return solveMosek(solver, pr);
}

struct Residuals {
    double stationarity{0};  // ||P x + q + A^T z||_inf
    double equality{0};      // ||(b - A x)_eq||_inf
    double primal_cone{0};   // max violation of s = b - A x in K
    double dual_cone{0};     // max violation of z in K*
    double complementarity{0};  // |s^T z|
    double objective{0};     // 1/2 x'Px + q'x
};

// Largest violation of v in K (nonneg part and second-order cones, zero part skipped).
inline double coneViolation(const Vec& v, const ConeDims& d) {
    double viol = 0;
    int i = d.zero;
    for (int k = 0; k < d.nonneg; ++k, ++i) viol = std::max(viol, -v[i]);
    for (int k : d.soc) {
        viol = std::max(viol, v.segment(i + 1, k - 1).norm() - v[i]);
        i += k;
    }
    return viol;
}

inline Residuals residuals(const Problem& pr, const Vec& x, const Vec& z) {
    Residuals r;
    const Vec s = pr.b - pr.A * x;
    const int p = pr.dims.zero;
    const int m = static_cast<int>(pr.A.rows()) - p;
    r.stationarity = (pr.P * x + pr.q + pr.A.transpose() * z).lpNorm<Eigen::Infinity>();
    r.equality = p ? s.head(p).lpNorm<Eigen::Infinity>() : 0.0;
    r.primal_cone = coneViolation(s, pr.dims);
    r.dual_cone = coneViolation(z, pr.dims);
    r.complementarity = m ? std::abs(s.tail(m).dot(z.tail(m))) : 0.0;
    r.objective = 0.5 * x.dot(pr.P * x) + pr.q.dot(x);
    return r;
}

/**
 * @brief Synthetic problem with the structure ClarabelAssembler produces for a Cardillo scene:
 * x = [v (6 per body); lambda_spring], P = diag(M, C/(theta dt^2)), equality rows
 * [Wg  -C] (springs between body pairs), then 3-row frictional contacts (normal row scaled by 1/mu)
 * between body pairs or body/ground; b is chosen so that the problem is feasible.
 */
inline Problem cardilloLikeProblem(int nbodies, int ncontacts, int nsprings, unsigned seed, int nfrictionless = 0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    std::uniform_int_distribution<int> body(0, nbodies - 1);
    const int nv = 6 * nbodies;
    const int n = nv + nsprings;
    const double dt = 1e-3, theta = 0.5;

    Problem pr;
    pr.dims.zero = nsprings;
    pr.dims.nonneg = nfrictionless;
    pr.dims.soc.assign(static_cast<std::size_t>(ncontacts), 3);
    const int m = nfrictionless + 3 * ncontacts;

    std::vector<Eigen::Triplet<double>> Pt, At;
    Vec v_old(nv);
    for (int i = 0; i < nv; ++i) {
        const double mass = (i % 6 < 3) ? 1.0 + 0.5 * (U(rng) + 1.0) : 0.01 + 0.005 * (U(rng) + 1.0);
        Pt.emplace_back(i, i, mass);
        v_old[i] = U(rng);
    }
    for (int k = 0; k < nsprings; ++k) Pt.emplace_back(nv + k, nv + k, (1e-6 * (1.5 + U(rng))) / (theta * dt * dt));

    auto bodyRow = [&](int row, int bod, double sign) {
        for (int c = 0; c < 6; ++c) At.emplace_back(row, 6 * bod + c, sign * U(rng));
    };
    for (int k = 0; k < nsprings; ++k) {
        const int b1 = body(rng);
        int b2 = body(rng);
        if (b2 == b1) b2 = (b1 + 1) % nbodies;
        bodyRow(k, b1, 1.0);
        bodyRow(k, b2, -1.0);
        At.emplace_back(k, nv + k, -Pt[nv + k].value());  // -C/(theta dt^2)
    }

    int row = nsprings;
    auto contactRows = [&](int nrows, const std::vector<double>& smu) {
        const int b1 = body(rng);
        const bool ground = (U(rng) > 0.3) || nbodies == 1;
        int b2 = body(rng);
        if (b2 == b1) b2 = (b1 + 1) % nbodies;
        for (int r = 0; r < nrows; ++r) {
            // -Smu W
            for (int c = 0; c < 6; ++c) At.emplace_back(row + r, 6 * b1 + c, -smu[r] * U(rng));
            if (!ground)
                for (int c = 0; c < 6; ++c) At.emplace_back(row + r, 6 * b2 + c, -smu[r] * U(rng));
        }
        row += nrows;
    };
    for (int k = 0; k < nfrictionless; ++k) contactRows(1, {1.0});
    for (int k = 0; k < ncontacts; ++k) {
        const double mu = 0.2 + 0.3 * (U(rng) + 1.0);
        contactRows(3, {1.0 / mu, 1.0, 1.0});
    }
    pr.P.resize(n, n);
    pr.P.setFromTriplets(Pt.begin(), Pt.end());
    pr.P.makeCompressed();
    pr.A.resize(nsprings + m, n);
    pr.A.setFromTriplets(At.begin(), At.end());
    pr.A.makeCompressed();

    // b = A x0 + s0 with s0 strictly inside the cones, so the problem is feasible (like a real
    // contact step) while most cones are still active at the optimum (x0 is far from it).
    Vec x0 = 0.1 * Vec::NullaryExpr(n, [&]() { return U(rng); });
    Vec s0(nsprings + m);
    s0.head(nsprings).setZero();
    for (int k = 0; k < nfrictionless; ++k) s0[nsprings + k] = 0.01 * (1.5 + U(rng));
    for (int k = 0, r = nsprings + nfrictionless; k < ncontacts; ++k, r += 3) {
        s0[r + 1] = 0.05 * U(rng);
        s0[r + 2] = 0.05 * U(rng);
        s0[r] = std::hypot(s0[r + 1], s0[r + 2]) + 0.01 * (1.5 + U(rng));
    }
    pr.b = pr.A * x0 + s0;
    pr.q = Vec::Zero(n);
    pr.q.head(nv) = -(pr.P.topLeftCorner(nv, nv) * v_old + dt * 9.81 * Vec::Ones(nv));
    return pr;
}

}  // namespace conic_test
