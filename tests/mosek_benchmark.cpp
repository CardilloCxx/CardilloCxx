// Benchmark of the MOSEK backend on synthetic Cardillo-shaped contact problems (see
// conic_test::cardilloLikeProblem), compared with Clarabel and QOCO on the same data.
//
// Reported per MOSEK configuration: model construction (setup), numerical update, optimize,
// solution extraction, total -- for the first solve and averaged over subsequent solves of a
// value-perturbed problem with an unchanged structure (the persistent-task path). Clarabel and QOCO
// rebuild their solver every step in Cardillo, so only their total time per solve is reported.
//
// Usage: mosek_benchmark [nbodies ncontacts nsprings repeats]

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <thread>

#include "conic_test_utils.hpp"

using namespace conic_test;
namespace mosek = cardillo::solver::mosek;

namespace {

struct Phases {
    double setup{0}, update{0}, optimize{0}, extract{0};
    int iters{0};
    double total() const { return setup + update + optimize + extract; }
};

Phases timedSolve(mosek::ConicSolver& s, const Problem& pr) {
    Phases ph;
    double t = now();
    if (s.structureMatches(pr.P, pr.A, pr.dims)) {
        s.update(pr.P, pr.q, pr.A, pr.b);
        ph.update = now() - t;
    } else {
        s.setup(pr.P, pr.q, pr.A, pr.b, pr.dims);
        ph.setup = now() - t;
    }
    t = now();
    const auto st = s.optimize();
    ph.optimize = now() - t;
    if (!mosek::isAcceptable(st)) std::cerr << "MOSEK: " << mosek::toString(st) << "\n";
    t = now();
    s.extract();
    ph.extract = now() - t;
    ph.iters = s.iterations();
    return ph;
}

void perturb(Problem& pr, std::mt19937& rng) {
    std::uniform_real_distribution<double> U(0.95, 1.05);
    for (int k = 0; k < pr.A.nonZeros(); ++k) pr.A.valuePtr()[k] *= U(rng);
    for (int k = 0; k < pr.b.size(); ++k) pr.b[k] *= U(rng);
    for (int k = 0; k < pr.q.size(); ++k) pr.q[k] *= U(rng);
}

void row(const char* label, const Phases& p) {
    std::printf("  %-22s setup %8.3f  update %8.3f  optimize %8.3f  extract %7.3f  total %8.3f ms  (%d it)\n", label, 1e3 * p.setup, 1e3 * p.update, 1e3 * p.optimize, 1e3 * p.extract,
                1e3 * p.total(), p.iters);
}

}  // namespace

int main(int argc, char** argv) {
    const int nbodies = argc > 1 ? std::atoi(argv[1]) : 2000;
    const int ncontacts = argc > 2 ? std::atoi(argv[2]) : 5000;
    const int nsprings = argc > 3 ? std::atoi(argv[3]) : 200;
    const int repeats = argc > 4 ? std::atoi(argv[4]) : 10;

    Problem base = cardilloLikeProblem(nbodies, ncontacts, nsprings, 1234, ncontacts / 20);
    std::printf("Cardillo-shaped problem: %d bodies, %d frictional + %d frictionless contacts, %d springs\n", nbodies, ncontacts, ncontacts / 20, nsprings);
    std::printf("  n = %ld variables, %ld constraint rows, nnz(A) = %ld\n\n", (long)base.P.cols(), (long)base.A.rows(), (long)base.A.nonZeros());

    // Reference solvers (full rebuild per solve, as in ClarabelSolver/QocoSolver).
    double tc = 0, tq = 0;
    Result rc, rq;
    for (int k = 0; k < 3; ++k) {
        rc = solveClarabel(base);
        tc += rc.seconds;
        rq = solveQoco(base);
        tq += rq.seconds;
    }
    std::printf("Clarabel: %8.3f ms per solve (setup + solve)%s\n", 1e3 * tc / 3, rc.ok ? "" : "  [NOT SOLVED]");
    std::printf("QOCO:     %8.3f ms per solve (setup + solve)%s\n\n", 1e3 * tq / 3, rq.ok ? "" : "  [NOT SOLVED]");

    std::vector<int> threads = {1, 2, 4, 8};
    const int hw = static_cast<int>(std::thread::hardware_concurrency());
    if (hw > 8) threads.push_back(hw);
    threads.push_back(0);  // MOSEK's own choice

    for (int nt : threads) {
        mosek::Options o;
        o.num_threads = nt;
        o.equality_duals = false;
        mosek::ConicSolver solver(o);
        Problem pr = base;
        std::mt19937 rng(99);

        const Phases first = timedSolve(solver, pr);
        Phases avg;
        for (int k = 0; k < repeats; ++k) {
            perturb(pr, rng);
            const Phases p = timedSolve(solver, pr);
            avg.setup += p.setup / repeats;
            avg.update += p.update / repeats;
            avg.optimize += p.optimize / repeats;
            avg.extract += p.extract / repeats;
            avg.iters += p.iters;
        }
        avg.iters /= std::max(1, repeats);

        // A fresh setup each step (what a changed contact set costs).
        Phases rebuild;
        for (int k = 0; k < repeats; ++k) {
            mosek::ConicSolver fresh(o);
            const Phases p = timedSolve(fresh, pr);
            rebuild.setup += p.setup / repeats;
            rebuild.optimize += p.optimize / repeats;
            rebuild.extract += p.extract / repeats;
            rebuild.iters += p.iters;
        }
        rebuild.iters /= std::max(1, repeats);

        std::printf("MOSEK, threads = %s\n", nt == 0 ? "auto" : std::to_string(nt).c_str());
        row("first solve", first);
        row("subsequent (update)", avg);
        row("subsequent (rebuild)", rebuild);

        if (nt == threads.front()) {
            const Result rm = solveMosek(solver, base);
            if (rm.ok && rc.ok) std::printf("  ||x_mosek - x_clarabel||_inf = %.2e, |f_mosek - f_clarabel| = %.2e\n", (rm.x - rc.x).lpNorm<Eigen::Infinity>(), std::abs(rm.obj - rc.obj));
        }
    }
    return 0;
}
