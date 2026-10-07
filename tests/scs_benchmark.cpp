// Benchmark of the SCS backend (direct, AMD + QDLDL) on synthetic Cardillo-shaped contact problems
// (see conic_test::cardilloLikeProblem), compared with Clarabel, QOCO, ConicXX and (when built)
// MOSEK on exactly the same data. The interior-point solvers run at tight tolerances (1e-9) and
// Clarabel's solution is the accuracy reference; SCS runs at eps in {1e-4, 1e-6}, cold and warm
// started. "warm" solves a sequence of value-perturbed problems (same structure, +-1%), each warm
// started from the previous solution, i.e. the time-stepping situation with a persistent contact set.
//
// Accuracy columns: rel. KKT = max(stationarity, primal/dual cone violation, complementarity) /
// max(1, ||q||, ||b||); dx = ||x - x_clarabel||_inf / max(1, ||x_clarabel||_inf).
//
// Usage: scs_benchmark [nbodies ncontacts nsprings repeats]

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "conic_test_utils.hpp"

using namespace conic_test;
namespace scs = cardillo::solver::scs;

namespace {

void row(const char* label, double ms, double iters, double kkt, double dx, const char* extra = "") {
    std::printf("  %-26s %10.3f ms  %8.1f it   rel.KKT %8.1e   dx %8.1e  %s\n", label, ms, iters, kkt, dx, extra);
}

void perturb(Problem& pr, std::mt19937& rng) {
    std::uniform_real_distribution<double> U(0.99, 1.01);
    for (int k = 0; k < pr.b.size(); ++k) pr.b[k] *= U(rng);
    for (int k = 0; k < pr.q.size(); ++k) pr.q[k] *= U(rng);
}

}  // namespace

int main(int argc, char** argv) {
    const int nbodies = argc > 1 ? std::atoi(argv[1]) : 500;
    const int ncontacts = argc > 2 ? std::atoi(argv[2]) : 1500;
    const int nsprings = argc > 3 ? std::atoi(argv[3]) : 50;
    const int repeats = argc > 4 ? std::atoi(argv[4]) : 10;

    const Problem base = cardilloLikeProblem(nbodies, ncontacts, nsprings, 1234, ncontacts / 20);
    std::printf("Cardillo-shaped problem: %d bodies, %d frictional + %d frictionless contacts, %d springs\n", nbodies, ncontacts, ncontacts / 20, nsprings);
    std::printf("  n = %ld variables, %ld constraint rows, nnz(A) = %ld\n\n", (long)base.P.cols(), (long)base.A.rows(), (long)base.A.nonZeros());

    // Interior-point references (median-free: best of 3 to suppress first-touch noise).
    auto best = [](auto&& fn) {
        Result r = fn();
        for (int k = 0; k < 2; ++k) {
            Result t = fn();
            if (t.seconds < r.seconds) r = t;
        }
        return r;
    };
    const Result rc = best([&] { return solveClarabel(base); });
    const Result rq = best([&] { return solveQoco(base); });
    const Result rx = best([&] { return solveConicxx(base); });
    std::printf("Interior point (tol 1e-9, setup + solve):\n");
    row("Clarabel", 1e3 * rc.seconds, std::nan(""), relKkt(base, rc), 0.0, rc.ok ? "" : "[NOT SOLVED]");
    row("QOCO", 1e3 * rq.seconds, std::nan(""), relKkt(base, rq), relDx(rq, rc), rq.ok ? "" : "[NOT SOLVED]");
    row("ConicXX", 1e3 * rx.seconds, std::nan(""), relKkt(base, rx), relDx(rx, rc), rx.ok ? "" : "[NOT SOLVED]");
#ifdef CARDILLO_HAVE_MOSEK
    const Result rm = best([&] { return solveMosek(base); });
    row("MOSEK", 1e3 * rm.seconds, std::nan(""), relKkt(base, rm), relDx(rm, rc), rm.ok ? "" : "[NOT SOLVED]");
#endif

    for (double eps : {1e-4, 1e-6}) {
        scs::Settings st;
        st.eps_abs = st.eps_rel = eps;
        std::printf("\nSCS direct, eps = %.0e (setup = scs_init incl. factorization):\n", eps);

        scs::Info info;
        const Result r = best([&] { return solveScs(base, st, nullptr, &info); });
        char extra[160];
        std::snprintf(extra, sizeof extra, "[%s] setup %.2f + solve %.2f ms (linsys %.2f, cone %.2f, AA %.2f), %d scale upd.", scs::toString(info.status), info.setup_ms, info.solve_ms,
                      info.lin_sys_ms, info.cone_ms, info.accel_ms, info.scale_updates);
        row("cold", 1e3 * r.seconds, info.iterations, relKkt(base, r), relDx(r, rc), extra);

        // Warm-started sequence vs cold-started sequence on the same perturbed problems.
        for (bool warm : {false, true}) {
            Problem pr = base;
            std::mt19937 rng(99);
            Result prev = r;
            double t = 0, it = 0, kkt = 0;
            scs::Settings stw = st;
            for (int k = 0; k < repeats; ++k) {
                perturb(pr, rng);
                scs::Info in;
                const Result rk = solveScs(pr, stw, warm ? &prev : nullptr, &in);
                t += rk.seconds / repeats;
                it += double(in.iterations) / repeats;
                kkt = std::max(kkt, relKkt(pr, rk));
                if (warm && in.scale > 0) stw.scale = in.scale;  // carry the adapted scale, as ScsSolver does
                prev = rk;
            }
            row(warm ? "sequence, warm (+scale)" : "sequence, cold", 1e3 * t, it, kkt, std::nan(""));
        }
    }
    return 0;
}
