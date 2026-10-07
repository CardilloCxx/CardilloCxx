// Tests of the MOSEK conic backend (cardillo::solver::mosek::ConicSolver) against analytic solutions
// and against Clarabel/QOCO on exactly the same problem data. Exit code = number of failed checks.

#include <cstdio>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

#include "conic_test_utils.hpp"

using namespace conic_test;
namespace mosek = cardillo::solver::mosek;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool cond, const std::string& what) {
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::cout << "    FAIL: " << what << "\n";
    }
}

void checkNear(double a, double b, double tol, const std::string& what) {
    const bool ok = std::abs(a - b) <= tol * std::max(1.0, std::abs(b));
    check(ok, what + " (" + std::to_string(a) + " vs " + std::to_string(b) + ")");
}

void checkVec(const Vec& a, const Vec& b, double tol, const std::string& what) {
    if (a.size() != b.size()) {
        check(false, what + " (size mismatch)");
        return;
    }
    const double err = (a - b).lpNorm<Eigen::Infinity>();
    const double scale = std::max(1.0, b.lpNorm<Eigen::Infinity>());
    check(err <= tol * scale, what + " (||diff||_inf = " + std::to_string(err) + ")");
}

void checkKkt(const Problem& pr, const Result& r, double tol, const std::string& who) {
    const Residuals res = residuals(pr, r.x, r.z);
    const double scale = std::max({1.0, pr.q.lpNorm<Eigen::Infinity>(), pr.b.lpNorm<Eigen::Infinity>()});
    check(res.stationarity <= tol * scale, who + " stationarity " + std::to_string(res.stationarity));
    check(res.equality <= tol * scale, who + " equality residual " + std::to_string(res.equality));
    check(res.primal_cone <= tol * scale, who + " primal cone violation " + std::to_string(res.primal_cone));
    check(res.dual_cone <= tol * scale, who + " dual cone violation " + std::to_string(res.dual_cone));
    check(res.complementarity <= tol * scale, who + " complementarity " + std::to_string(res.complementarity));
    checkNear(r.obj, res.objective, tol, who + " reported objective vs 1/2x'Px+q'x");
}

// Interior-point solutions stopped at a relative gap of ~1e-8 (MOSEK's default) are accurate to
// roughly sqrt(gap) in x and z, and the references disagree among themselves at that level, so
// vectors are compared across solvers with this tolerance; KKT residuals and objectives (which
// converge like the gap) are checked with the tight per-test tolerance.
constexpr double kCrossSolverVecTol = 1e-4;

/// Solves with MOSEK, Clarabel and QOCO; checks MOSEK's KKT residuals and agreement of x/obj
/// (and of the duals when `compare_duals`, i.e. when they are unique).
Result crossCheck(const Problem& pr, double tol, bool compare_duals = true) {
    const Result rm = solveMosek(pr);
    const Result rc = solveClarabel(pr);
    const Result rq = solveQoco(pr);
    check(rm.ok, "MOSEK solved");
    check(rc.ok, "Clarabel solved");
    check(rq.ok, "QOCO solved");
    if (!rm.ok) return rm;
    checkKkt(pr, rm, tol, "MOSEK");
    if (rc.ok) {
        checkVec(rm.x, rc.x, kCrossSolverVecTol, "x MOSEK vs Clarabel");
        checkNear(rm.obj, rc.obj, tol, "objective MOSEK vs Clarabel");
        if (compare_duals) checkVec(rm.z, rc.z, kCrossSolverVecTol, "z MOSEK vs Clarabel");
    }
    if (rq.ok) {
        checkVec(rm.x, rq.x, kCrossSolverVecTol, "x MOSEK vs QOCO");
        checkNear(rm.obj, rq.obj, tol, "objective MOSEK vs QOCO");
    }
    return rm;
}

SpMat sparse(const Eigen::MatrixXd& D) {
    SpMat S = D.sparseView();
    S.makeCompressed();
    return S;
}

SpMat identity(int n) {
    SpMat I(n, n);
    I.setIdentity();
    I.makeCompressed();
    return I;
}

// ---------------------------------------------------------------------------------------------

void testUnconstrainedDiagonal() {
    // min 1/2 (2 x0^2 + 4 x1^2) - 2 x0 - 8 x1  ->  x = (1, 2), f = -1 - 8 = -9
    Problem pr;
    pr.P = sparse(Eigen::Vector2d(2, 4).asDiagonal().toDenseMatrix());
    pr.q = Vec::Zero(2);
    pr.q << -2, -8;
    pr.A.resize(0, 2);
    pr.b.resize(0);
    const Result r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector2d(1, 2), 1e-6, "analytic x");
    checkNear(r.obj, -9.0, 1e-6, "analytic objective (factor 1/2 convention)");
}

void testUnconstrainedGeneralP() {
    // Non-diagonal P = [[4,1],[1,2]], q = (1,1): x = -P^{-1} q = -(1,3)/7.
    Eigen::Matrix2d Pd;
    Pd << 4, 1, 1, 2;
    Problem pr;
    pr.P = sparse(Pd);
    pr.q = Vec::Ones(2);
    pr.A.resize(0, 2);
    pr.b.resize(0);
    const Eigen::Vector2d xs(-1.0 / 7, -3.0 / 7);
    const Result r = crossCheck(pr, 1e-6);
    checkVec(r.x, xs, 1e-6, "analytic x (full symmetric P)");
    checkNear(r.obj, 0.5 * xs.dot(Pd * xs) + pr.q.dot(xs), 1e-6, "analytic objective");

    // Only the lower triangle must be read: passing just the lower triangle gives the same result.
    Problem lower = pr;
    lower.P = pr.P.triangularView<Eigen::Lower>();
    lower.P.makeCompressed();
    const Result rl = solveMosek(lower);
    check(rl.ok, "MOSEK solved (lower-triangular P)");
    if (rl.ok) checkVec(rl.x, xs, 1e-6, "analytic x (lower-triangular P)");

    // Larger random SPD P exercises the LDL^T permutation.
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(-1, 1);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(8, 8);
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            if (i == j || U(rng) > 0.5) B(i, j) = U(rng);
    Eigen::MatrixXd P8 = B * B.transpose() + 0.1 * Eigen::MatrixXd::Identity(8, 8);
    Problem big;
    big.P = sparse(P8);
    big.q = Vec::NullaryExpr(8, [&]() { return U(rng); });
    big.A = -identity(8);  // x >= 0 so the cones are active as well
    big.b = Vec::Zero(8);
    big.dims.nonneg = 8;
    crossCheck(big, 1e-6);
}

void testEqualityQp() {
    // min 1/2 ||x||^2  s.t. x0 + x1 + x2 = 3  ->  x = (1,1,1), z_eq = -1
    Problem pr;
    pr.P = identity(3);
    pr.q = Vec::Zero(3);
    pr.A = sparse(Eigen::RowVector3d(1, 1, 1));
    pr.b = Vec::Constant(1, 3.0);
    pr.dims.zero = 1;
    const Result r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector3d(1, 1, 1), 1e-6, "analytic x");
    checkVec(r.z, Vec::Constant(1, -1.0), 1e-6, "analytic equality dual");
    checkNear(r.obj, 1.5, 1e-6, "analytic objective");
}

void testNonnegativeCone() {
    // min 1/2 ||x - (1,-2)||^2  s.t. x >= 0  (A = -I, b = 0)  ->  x = (1, 0), z = (0, 2)
    Problem pr;
    pr.P = identity(2);
    pr.q = Vec::Zero(2);
    pr.q << -1, 2;
    pr.A = -identity(2);
    pr.b = Vec::Zero(2);
    pr.dims.nonneg = 2;
    const Result r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector2d(1, 0), 1e-6, "analytic x");
    checkVec(r.z, Eigen::Vector2d(0, 2), 1e-6, "analytic z");
}

void testSingleSoc() {
    // Projection of a = (0,3,4) onto Q^3: x = ((0+5)/2)(1, 3/5, 4/5) = (2.5, 1.5, 2).
    Problem pr;
    pr.P = identity(3);
    pr.q = -Eigen::Vector3d(0, 3, 4);
    pr.A = -identity(3);
    pr.b = Vec::Zero(3);
    pr.dims.soc = {3};
    Result r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector3d(2.5, 1.5, 2.0), 1e-6, "analytic x (cone x in Q)");
    // Dual z = x - a (stationarity x - a - z = 0), lies on the boundary of Q.
    checkVec(r.z, Eigen::Vector3d(2.5, -1.5, -2.0), 1e-6, "analytic z");

    // Sign test for g = h: b = -(1,0,0) gives s = x - e0 in Q, i.e. x0 - 1 >= ||(x1,x2)||.
    // Projection of a - e0 = (-1,3,4): ((-1+5)/2)(1, .6, .8) = (2, 1.2, 1.6)  ->  x = (3, 1.2, 1.6).
    pr.b = -Eigen::Vector3d(1, 0, 0);
    r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector3d(3.0, 1.2, 1.6), 1e-6, "analytic x (shifted cone, sign of h)");

    // Sign test for F = -G: G = +I (s = -x in Q) projects a onto -Q: x = -proj_Q(-a) = (-2.5, 1.5, 2)
    // only if the solver honours  h - Gx in K  (with F = +G it would give (2.5, 1.5, 2) again).
    pr.A = identity(3);
    pr.b = Vec::Zero(3);
    pr.q = -Eigen::Vector3d(0, 3, 4);
    r = crossCheck(pr, 1e-6);
    checkVec(r.x, Eigen::Vector3d(-2.5, 1.5, 2.0), 1e-6, "analytic x (G = +I, sign of F)");
}

void testMultipleCones() {
    // x in R^9; rows: 1 nonneg (x8 >= 0.5), Q^3 on (x0,x1,x2), Q^4 on (x6,x3,x4,x5) (permuted
    // variables), Q^3 on (x7, x0, x3). Targets chosen so that all cones are active.
    const int n = 9;
    Problem pr;
    pr.P = identity(n);
    Vec a(n);
    a << 0, 3, 4, 1, -2, 2, 0.5, -1, -3;
    pr.q = -a;
    std::vector<Eigen::Triplet<double>> t;
    int row = 0;
    t.emplace_back(row++, 8, -1.0);  // s = -0.5 + x8 >= 0  (b = -0.5)
    for (int v : {0, 1, 2}) t.emplace_back(row++, v, -1.0);
    for (int v : {6, 3, 4, 5}) t.emplace_back(row++, v, -1.0);
    for (int v : {7, 0, 3}) t.emplace_back(row++, v, -1.0);
    pr.A.resize(row, n);
    pr.A.setFromTriplets(t.begin(), t.end());
    pr.A.makeCompressed();
    pr.b = Vec::Zero(row);
    pr.b[0] = -0.5;
    pr.dims.nonneg = 1;
    pr.dims.soc = {3, 4, 3};
    const Result r = crossCheck(pr, 1e-6);
    if (r.ok) {
        check(r.x[8] >= 0.5 - 1e-7, "nonneg cone respected");
        const Vec s = pr.b - pr.A * r.x;
        check(std::abs(s[1] - s.segment(2, 2).norm()) < 1e-6, "cone 1 active");
        check(std::abs(s[4] - s.segment(5, 3).norm()) < 1e-6, "cone 2 (dim 4) active");
    }
}

void testEqualityAndSoc() {
    // Small Cardillo-shaped problem: springs (equalities) + frictionless + frictional contacts.
    for (unsigned seed : {1u, 2u, 3u}) {
        const Problem pr = cardilloLikeProblem(4, 5, 2, seed, 2);
        crossCheck(pr, 1e-5);
    }
}

void testLargeSparse() {
    // Same structure as a Cardillo contact scene at a scale of a few hundred bodies.
    const Problem pr = cardilloLikeProblem(400, 900, 60, 42, 30);
    std::cout << "    n = " << pr.P.cols() << ", rows = " << pr.A.rows() << ", nnz(A) = " << pr.A.nonZeros() << "\n";
    const Result r = crossCheck(pr, 1e-5, /*compare_duals=*/false);
    (void)r;
}

void testEpigraphLayouts() {
    // The single-cone and separable epigraphs of 1/2 x'Px model the same problem.
    const Problem pr = cardilloLikeProblem(10, 12, 2, 21, 2);
    mosek::Options single, separable;
    single.quad_model = mosek::Options::QuadModel::Single;
    separable.quad_model = mosek::Options::QuadModel::Separable;
    mosek::ConicSolver s1(single), s2(separable);
    const Result r1 = solveMosek(s1, pr), r2 = solveMosek(s2, pr);
    check(r1.ok && r2.ok, "both layouts solve");
    if (r1.ok && r2.ok) {
        checkVec(r1.x, r2.x, kCrossSolverVecTol, "x single vs separable");
        checkNear(r1.obj, r2.obj, 1e-6, "objective single vs separable");
        checkKkt(pr, r1, 1e-5, "MOSEK (single epigraph)");
    }
    // Without cone rows Auto passes P natively (MSK_putqobj); the epigraph models must agree.
    Problem eq = cardilloLikeProblem(10, 0, 8, 22, 0);
    mosek::Options native;
    native.quad_model = mosek::Options::QuadModel::Native;
    for (const auto& o : {native, single, separable}) {
        mosek::ConicSolver s(o);
        const Result r = solveMosek(s, eq);
        check(r.ok, "equality-only problem solves");
        if (r.ok) checkKkt(eq, r, 1e-6, "MOSEK (equality-only)");
    }
    mosek::ConicSolver bad(native);
    bool threw = false;
    try {
        bad.setup(pr.P, pr.q, pr.A, pr.b, pr.dims);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "native quadratic objective with cone rows is rejected");

    // Switching the option invalidates the cached structure.
    s1.setOptions(separable);
    check(!s1.structureMatches(pr.P, pr.A, pr.dims), "epigraph option change forces a rebuild");
}

void testScaling() {
    // Badly scaled Cardillo-shaped problem (inertias/compliances over many decades): the internal
    // pre-scaling must be invisible to the caller, i.e. x, z and the objective are unchanged.
    Problem pr = cardilloLikeProblem(20, 25, 6, 31, 3);
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> E(-8, 6);
    Vec D(pr.P.cols());
    for (int j = 0; j < D.size(); ++j) D[j] = std::pow(10.0, E(rng) / 2);
    pr.P = D.asDiagonal() * pr.P * D.asDiagonal();  // same problem in x / D, now poorly scaled
    pr.q = D.cwiseProduct(pr.q);
    pr.A = pr.A * D.asDiagonal();
    pr.P.makeCompressed();
    pr.A.makeCompressed();
    mosek::Options off;
    off.scaling = false;
    mosek::ConicSolver s_on, s_off(off);
    const Result on = solveMosek(s_on, pr);
    check(on.ok, "scaled solve");
    if (on.ok) checkKkt(pr, on, 1e-5, "MOSEK (scaled)");
    const Result un = solveMosek(s_off, pr);
    if (on.ok && un.ok) checkNear(on.obj, un.obj, 1e-6, "objective scaled vs unscaled");
    const Result rc = solveClarabel(pr);
    if (on.ok && rc.ok) {
        checkNear(on.obj, rc.obj, 1e-6, "objective scaled MOSEK vs Clarabel");
        checkVec(on.z, rc.z, kCrossSolverVecTol, "z scaled MOSEK vs Clarabel");
    }
}

void testUpdatePath() {
    // Same structure, new values: the update path must give the same result as a fresh setup.
    Problem pr = cardilloLikeProblem(30, 40, 5, 11, 3);
    mosek::ConicSolver persistent;
    check(!persistent.structureMatches(pr.P, pr.A, pr.dims), "no structure before first solve");
    Result r0 = solveMosek(persistent, pr);
    check(r0.ok, "first solve");
    check(persistent.structureMatches(pr.P, pr.A, pr.dims), "structure cached after first solve");

    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(0.8, 1.2);
    for (int step = 0; step < 3; ++step) {
        for (int k = 0; k < pr.A.nonZeros(); ++k) pr.A.valuePtr()[k] *= U(rng);
        for (int k = 0; k < pr.P.nonZeros(); ++k) pr.P.valuePtr()[k] *= U(rng);
        pr.q *= U(rng);
        pr.b = pr.b.cwiseProduct(Vec::NullaryExpr(pr.b.size(), [&]() { return U(rng); }));
        check(persistent.structureMatches(pr.P, pr.A, pr.dims), "structure still matches after value change");
        const bool rebuilt = persistent.load(pr.P, pr.q, pr.A, pr.b, pr.dims);
        check(!rebuilt, "value change uses update path");
        const Result ru = solveMosek(persistent, pr);
        const Result rf = solveMosek(pr);
        check(ru.ok && rf.ok, "update/fresh solves succeed");
        if (ru.ok && rf.ok) {
            checkVec(ru.x, rf.x, 1e-7, "update path x == fresh setup x");
            checkVec(ru.z, rf.z, 1e-6, "update path z == fresh setup z");
        }
        checkKkt(pr, ru, 1e-5, "MOSEK (update path)");
    }

    // A changed structure must trigger a rebuild and still be correct.
    const Problem pr2 = cardilloLikeProblem(30, 41, 5, 12, 3);
    check(!persistent.structureMatches(pr2.P, pr2.A, pr2.dims), "changed structure detected");
    const bool rebuilt = persistent.load(pr2.P, pr2.q, pr2.A, pr2.b, pr2.dims);
    check(rebuilt, "changed structure rebuilds the task");
    const Result r2 = solveMosek(persistent, pr2);
    const Result r2c = solveClarabel(pr2);
    if (r2.ok && r2c.ok) checkVec(r2.x, r2c.x, kCrossSolverVecTol, "rebuilt task x vs Clarabel");

    // Changing only the cone composition (same row count and pattern) must also rebuild.
    Problem pr3 = pr2;
    pr3.dims.nonneg += 3;
    pr3.dims.soc.pop_back();
    check(!persistent.structureMatches(pr3.P, pr3.A, pr3.dims), "changed cone dims detected");
}

void testStatusAndErrors() {
    // Primal infeasible: x >= 1 and -x >= 1.
    Problem pr;
    pr.P = identity(1);
    pr.q = Vec::Zero(1);
    pr.A = sparse((Eigen::MatrixXd(2, 1) << -1, 1).finished());
    pr.b = Vec::Constant(2, -1.0);
    pr.dims.nonneg = 2;
    mosek::ConicSolver s;
    const auto st = s.solve(pr.P, pr.q, pr.A, pr.b, pr.dims);
    check(st == mosek::Status::PrimalInfeasible, std::string("infeasible problem reported as ") + mosek::toString(st));

    // Dual infeasible (unbounded): min -x s.t. x >= 0, no quadratic term.
    Problem unb;
    unb.P.resize(1, 1);
    unb.q = Vec::Constant(1, -1.0);
    unb.A = -identity(1);
    unb.b = Vec::Zero(1);
    unb.dims.nonneg = 1;
    const auto su = s.solve(unb.P, unb.q, unb.A, unb.b, unb.dims);
    check(su == mosek::Status::DualInfeasible, std::string("unbounded problem reported as ") + mosek::toString(su));

    // Iteration limit.
    mosek::Options o;
    o.max_iterations = 1;
    mosek::ConicSolver limited(o);
    const Problem big = cardilloLikeProblem(20, 30, 2, 3, 0);
    const auto sl = limited.solve(big.P, big.q, big.A, big.b, big.dims);
    check(sl == mosek::Status::IterationLimit, std::string("iteration-limited solve reported as ") + mosek::toString(sl));

    auto throws = [](const std::function<void()>& f) {
        try {
            f();
        } catch (const std::invalid_argument&) {
            return true;
        }
        return false;
    };
    // Inconsistent cone dimensions / non-convex P are rejected.
    Problem bad = pr;
    bad.dims.nonneg = 1;
    check(throws([&] { s.setup(bad.P, bad.q, bad.A, bad.b, bad.dims); }), "inconsistent cone dims throw");
    Problem ncvx = pr;
    ncvx.P = -identity(1);
    check(throws([&] { s.setup(ncvx.P, ncvx.q, ncvx.A, ncvx.b, ncvx.dims); }), "negative P throws");

    // MOSEK API errors are converted into exceptions carrying MOSEK's code name.
    try {
        mosek::check(MSK_RES_ERR_INDEX, "MSK_test");
        check(false, "mosek::check throws");
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        check(msg.find("MSK_test") != std::string::npos && msg.find("MSK_RES_ERR_INDEX") != std::string::npos, "error message names operation and code: " + msg);
    }
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, void (*)()>> tests = {
        {"unconstrained QP (diagonal P)", testUnconstrainedDiagonal},
        {"unconstrained QP (general P)", testUnconstrainedGeneralP},
        {"equality constrained QP", testEqualityQp},
        {"nonnegative cone", testNonnegativeCone},
        {"single SOC (+ sign tests)", testSingleSoc},
        {"multiple cones", testMultipleCones},
        {"equality + SOC (Cardillo-shaped)", testEqualityAndSoc},
        {"large sparse (Cardillo-shaped)", testLargeSparse},
        {"epigraph layouts", testEpigraphLayouts},
        {"pre-scaling", testScaling},
        {"persistent task / update path", testUpdatePath},
        {"status mapping and errors", testStatusAndErrors},
    };
    for (const auto& [name, fn] : tests) {
        const int before = g_failures;
        std::cout << "[ RUN  ] " << name << std::endl;
        try {
            fn();
        } catch (const std::exception& e) {
            check(false, std::string("exception: ") + e.what());
        }
        std::cout << (g_failures == before ? "[  OK  ] " : "[ FAIL ] ") << name << std::endl;
    }
    std::cout << "\n" << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures;
}
