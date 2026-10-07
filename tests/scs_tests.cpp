// Tests of the SCS conic backend (cardillo::solver::scs::Solver, direct linear-system backend)
// against analytic solutions and against Clarabel on exactly the same problem data. Exit code =
// number of failed checks.

#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "conic_test_utils.hpp"

using namespace conic_test;
namespace scs = cardillo::solver::scs;

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

void checkVec(const Vec& a, const Vec& b, double tol, const std::string& what) {
    if (a.size() != b.size()) {
        check(false, what + " (size mismatch)");
        return;
    }
    const double err = (a - b).lpNorm<Eigen::Infinity>();
    const double scale = std::max(1.0, b.lpNorm<Eigen::Infinity>());
    check(err <= tol * scale, what + " (||diff||_inf = " + std::to_string(err) + ")");
}

scs::Settings tight() {
    scs::Settings s;
    s.eps_abs = 1e-9;
    s.eps_rel = 1e-9;
    return s;
}

SpMat sparse(const Eigen::MatrixXd& M) { return M.sparseView(); }

// x = (t, x1, x2): min 1/2 ||x||^2 - 2 x1  s.t.  t = 1 (zero cone), (t, x1, x2) in Q^3.
// With t = 1 the x1 part 1/2 x1^2 - 2 x1 is minimized at x1 = 2, which the cone caps at x1 = 1, so
// x* = (1, 1, 0).
void testTinySoc() {
    Problem pr;
    pr.P = sparse(Eigen::MatrixXd::Identity(3, 3));
    pr.q = Vec(3);
    pr.q << 0, -2, 0;  // pulls x1 to 2, the cone caps it at t = 1
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 3);
    A(0, 0) = 1;  // t = 1
    A.block(1, 0, 3, 3) = -Eigen::MatrixXd::Identity(3, 3);  // s = x in Q^3
    pr.A = sparse(A);
    pr.b = Vec::Zero(4);
    pr.b[0] = 1;
    pr.dims.zero = 1;
    pr.dims.soc = {3};

    const Result r = solveScs(pr, tight());
    check(r.ok, "SCS solved the tiny SOCP");
    Vec x_ref(3);
    x_ref << 1, 1, 0;
    checkVec(r.x, x_ref, 1e-6, "tiny SOCP x");
    // Stationarity P x + q + A'y = 0 with A = [e1'; -I]: y0 - y1 = -1, y2 = -1, y3 = 0; y_soc =
    // (y1, y2, y3) in Q^3 with s'y_soc = y1 - 1 = 0 (s = x* on the boundary) -> y = (0, 1, -1, 0).
    Vec y_ref(4);
    y_ref << 0, 1, -1, 0;
    checkVec(r.z, y_ref, 1e-6, "tiny SOCP y");
    const Residuals res = residuals(pr, r.x, r.z);
    check(res.stationarity < 1e-7 && res.primal_cone < 1e-7 && res.dual_cone < 1e-7 && res.complementarity < 1e-7, "tiny SOCP KKT residuals");
}

// Nonnegative cone: min 1/2||x - a||^2 s.t. x >= 0 -> x = max(a, 0).
void testNonneg() {
    Problem pr;
    const int n = 5;
    pr.P = sparse(Eigen::MatrixXd::Identity(n, n));
    Vec a(n);
    a << 1.0, -2.0, 0.5, -0.1, 3.0;
    pr.q = -a;
    pr.A = sparse(-Eigen::MatrixXd::Identity(n, n));
    pr.b = Vec::Zero(n);
    pr.dims.nonneg = n;
    const Result r = solveScs(pr, tight());
    check(r.ok, "SCS solved the nonnegative QP");
    checkVec(r.x, a.cwiseMax(0.0), 1e-7, "nonnegative QP x = max(a, 0)");
}

// Cardillo-shaped problems (springs, frictionless and frictional contacts): SCS vs Clarabel.
void testCardilloShaped(int nbodies, int ncontacts, int nsprings, unsigned seed) {
    const Problem pr = cardilloLikeProblem(nbodies, ncontacts, nsprings, seed, ncontacts / 10);
    const Result rc = solveClarabel(pr);
    check(rc.ok, "Clarabel reference solved");
    scs::Info info;
    const Result rs = solveScs(pr, tight(), nullptr, &info);
    check(rs.ok, std::string("SCS solved (") + scs::toString(info.status) + ", " + std::to_string(info.iterations) + " it)");
    if (!rs.ok || !rc.ok) return;
    checkVec(rs.x, rc.x, 1e-5, "x vs Clarabel");
    checkVec(rs.z, rc.z, 1e-4, "y vs Clarabel z");
    const double fscale = std::max(1.0, std::abs(rc.obj));
    check(std::abs(rs.obj - rc.obj) <= 1e-6 * fscale, "objective vs Clarabel (" + std::to_string(rs.obj) + " vs " + std::to_string(rc.obj) + ")");
}

void testCardilloSmall() { testCardilloShaped(10, 20, 3, 7); }
void testCardilloMedium() { testCardilloShaped(200, 500, 20, 11); }

// Warm starting from the exact solution must terminate almost immediately, and from the solution of
// a slightly perturbed problem must need fewer iterations than a cold start.
void testWarmStart() {
    const Problem pr = cardilloLikeProblem(100, 250, 10, 3, 10);
    scs::Settings st = tight();
    st.eps_abs = st.eps_rel = 1e-7;
    scs::Info cold, hot, perturbed;
    const Result r0 = solveScs(pr, st, nullptr, &cold);
    check(r0.ok, "cold solve");
    const Result r1 = solveScs(pr, st, &r0, &hot);
    check(r1.ok, "warm solve from the solution");
    check(hot.iterations <= 50, "warm start from the solution converges quickly (" + std::to_string(hot.iterations) + " vs cold " + std::to_string(cold.iterations) + ")");

    Problem pr2 = pr;
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(0.99, 1.01);
    for (int k = 0; k < pr2.b.size(); ++k) pr2.b[k] *= U(rng);
    for (int k = 0; k < pr2.q.size(); ++k) pr2.q[k] *= U(rng);
    scs::Info cold2;
    const Result r2c = solveScs(pr2, st, nullptr, &cold2);
    const Result r2w = solveScs(pr2, st, &r0, &perturbed);
    check(r2c.ok && r2w.ok, "perturbed problem solved cold and warm");
    check(perturbed.iterations < cold2.iterations, "warm start on a perturbed problem saves iterations (" + std::to_string(perturbed.iterations) + " vs " + std::to_string(cold2.iterations) + ")");
}

void testProjection() {
    scs::ConeDims d;
    d.zero = 1;
    d.nonneg = 2;
    d.soc = {3, 3, 3};
    Vec v(12);
    v << 5.0,                // zero cone: untouched
        -1.0, 2.0,           // nonneg
        2.0, 1.0, 0.0,       // inside
        -3.0, 1.0, 1.0,      // in the polar cone -> 0
        0.0, 3.0, 4.0;       // boundary projection: ((0 + 5)/2, (3,4) * 2.5/5)
    scs::projectOntoCones(v, d);
    Vec ref(12);
    ref << 5.0, 0.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.5, 1.5, 2.0;
    checkVec(v, ref, 1e-14, "cone projection");
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, void (*)()>> tests = {
        {"cone projection", testProjection},
        {"tiny SOCP (analytic)", testTinySoc},
        {"nonnegative QP (analytic)", testNonneg},
        {"Cardillo-shaped small vs Clarabel", testCardilloSmall},
        {"Cardillo-shaped medium vs Clarabel", testCardilloMedium},
        {"warm start", testWarmStart},
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
