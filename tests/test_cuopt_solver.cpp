// cuOpt backend (cuopt_backend.hpp) on small cone programs with closed-form or known solutions. Needs
// an NVIDIA GPU: without one the test exits with 77 (ctest: skipped), with CARDILLO_REQUIRE_GPU=1 it
// fails instead.
//  1. factor 1/2 of cuOpt's objective: min x^2 - 2x s.t. x <= 10 gives x = 1 (x = 1/2 if Q = P had
//     been passed);
//  2. projections onto L^3, L^2 and L^1: inside, onto the apex (polar cone, the case of a sticking
//     contact) and onto the boundary;
//  3. random problems with a known KKT point mixing zero, nonnegative and second-order cones, with
//     left-out variables like the percussions of perfect constraints: x, A'z, the optimal value and
//     the KKT residuals (z itself only where it is unique);
//  4. an iteration limit of 1 is not reported as optimal, an infeasible problem fails, and the next
//     solve is unaffected.

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "cuopt_gpu.hpp"
#include "cuopt_test_problems.hpp"

using namespace cuopt_test;
using namespace cardillo::solver::cuopt;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) {
        ++g_failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

constexpr real_t kTol = 1e-8;     // tolerances of the barrier method
constexpr real_t kAccept = 1e-6;  // accepted relative error of the solution and the KKT residuals: 100 times kTol

Settings testSettings() {
    Settings s;
    s.tol_abs = kTol;
    s.tol_rel = kTol;
    s.deterministic = true;
    return s;
}

real_t relError(const VectorXr& a, const VectorXr& ref) {
    return (a - ref).lpNorm<Eigen::Infinity>() / (1 + ref.lpNorm<Eigen::Infinity>());
}

void solveKnown(Solver& solver, const KnownProblem& pr, const std::string& tag, bool unique_z) {
    solver.setup(pr.P, pr.q, pr.A, pr.b, pr.dims);
    const Status status = solver.solve();
    const Info& info = solver.info();
    std::cout << tag << ": " << toString(status) << ", " << info.iterations << " iterations, " << info.solve_time * 1e3 << " ms" << std::endl;
    check(status == Status::Optimal, tag + ": optimal (" + info.message + ")");
    if (status == Status::Failed) return;

    solver.recover(pr.A, pr.b);
    const VectorXr& x = solver.x();
    const VectorXr& z = solver.z();
    check(relError(x, pr.x) < kAccept, tag + ": x (error " + std::to_string(relError(x, pr.x)) + ")");
    const VectorXr ATz = pr.A.transpose() * z;
    const VectorXr ATz_ref = pr.A.transpose() * pr.z;
    check(relError(ATz, ATz_ref) < kAccept, tag + ": A'z (error " + std::to_string(relError(ATz, ATz_ref)) + ")");
    if (unique_z) check(relError(z, pr.z) < kAccept, tag + ": z (error " + std::to_string(relError(z, pr.z)) + ")");
    check(std::abs(info.objective + pr.value) < kAccept * (1 + std::abs(pr.value)), tag + ": objective " + std::to_string(info.objective) + " = -optimal value " + std::to_string(-pr.value));
    const KktResiduals r = kktResiduals(pr.P, pr.q, pr.A, pr.b, pr.dims, x, solver.s(), z);
    check(r.max() < kAccept, tag + ": KKT residuals primal " + std::to_string(r.primal) + ", dual " + std::to_string(r.dual) + ", cone(z) " + std::to_string(r.dual_cone) + ", complementarity " +
                                 std::to_string(r.complementarity));
}

KnownProblem halfFactorProblem() {
    KnownProblem pr;
    pr.dims.nonneg = 1;
    pr.P_full.resize(1, 1);
    pr.P_full.insert(0, 0) = 2;
    pr.P = pr.P_full;
    pr.A.resize(1, 1);
    pr.A.insert(0, 0) = 1;
    pr.q = VectorXr::Constant(1, -2);
    pr.b = VectorXr::Constant(1, 10);
    pr.x = VectorXr::Constant(1, 1);
    pr.s = VectorXr::Constant(1, 9);
    pr.z = VectorXr::Zero(1);
    pr.value = -1;
    return pr;
}

// x <= -1 and x >= 1.
KnownProblem infeasibleProblem() {
    KnownProblem pr;
    pr.dims.nonneg = 2;
    pr.P_full.resize(1, 1);
    pr.P_full.insert(0, 0) = 1;
    pr.P = pr.P_full;
    pr.A.resize(2, 1);
    pr.A.insert(0, 0) = 1;
    pr.A.insert(1, 0) = -1;
    pr.q = VectorXr::Zero(1);
    pr.b = VectorXr::Constant(2, -1);
    return pr;
}

}  // namespace

int main() {
    if (!gpuAvailable("test_cuopt_solver")) return gpuRequired() ? 1 : 77;
    std::cout << "cuOpt " << Solver::version() << std::endl;

    try {
        Solver solver(testSettings());
        solveKnown(solver, halfFactorProblem(), "factor 1/2", true);

        for (const VectorXr& p : {VectorXr(Vector3r(2, 1, 0)), VectorXr(Vector3r(-2, 1, 0)), VectorXr(Vector3r(0, 1, 0)), VectorXr(Eigen::Vector2d(0, 1)), VectorXr(Eigen::Vector2d(-1, 3)),
                                  VectorXr(VectorXr::Constant(1, -1)), VectorXr(VectorXr::Constant(1, 2))}) {
            std::string tag = "projection of (";
            for (Eigen::Index i = 0; i < p.size(); ++i) tag += (i ? ", " : "") + std::to_string(p[i]);
            solveKnown(solver, coneProjection(p), tag + ")", true);
        }

        ConeDims mixed;
        mixed.zero = 3;
        mixed.nonneg = 4;
        mixed.soc = {3, 3, 3, 3, 2, 1, 3, 2, 1};
        solveKnown(solver, randomProblem(1, 12, 2, mixed, true), "mixed, diagonal P", false);
        solveKnown(solver, randomProblem(2, 12, 2, mixed, false), "mixed, coupled P", false);
        ConeDims contacts;
        contacts.soc.assign(4, 3);
        solveKnown(solver, randomProblem(3, 18, 2, contacts, true), "frictional contacts, unique z", true);

        Settings limited_settings = testSettings();
        limited_settings.iteration_limit = 1;
        Solver limited(limited_settings);
        const KnownProblem pr = randomProblem(2, 12, 2, mixed, false);
        limited.setup(pr.P, pr.q, pr.A, pr.b, pr.dims);
        const Status limited_status = limited.solve();
        check(limited_status != Status::Optimal, std::string("iteration limit 1 is not optimal (") + toString(limited_status) + ")");

        const KnownProblem infeasible = infeasibleProblem();
        solver.setup(infeasible.P, infeasible.q, infeasible.A, infeasible.b, infeasible.dims);
        const Status infeasible_status = solver.solve();
        std::cout << "infeasible problem: " << toString(infeasible_status) << ": " << solver.info().message << std::endl;
        check(infeasible_status == Status::Failed, "an infeasible problem fails");
        solveKnown(solver, coneProjection(Vector3r(0, 1, 0)), "projection after a failed solve", true);
    } catch (const std::exception& ex) {
        check(false, std::string("exception: ") + ex.what());
    }
    if (g_failures == 0) std::cout << "test_cuopt_solver: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
