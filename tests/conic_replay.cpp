// Replays cone programs dumped from a real scene by ScsSolver (scs.dump_dir) through every conic
// backend on exactly the same data, so that time and accuracy are compared problem by problem
// instead of along diverging trajectories:
//
//   * interior point: Clarabel, QOCO, ConicXX and (when built) MOSEK at tolerance `ipm_tol`;
//   * SCS direct at eps in {1e-4, 1e-6}, cold and warm started. "warm" uses exactly the warm start
//     (x0, y0, s0, initial scale) ScsSolver built for that step from the tracked contact impulses.
//
// Times are wall clock per problem including setup (all backends set up from scratch every step in
// Cardillo, apart from MOSEK/ConicXX's update path, which is not exercised here). The reference for
// dx is Clarabel at 1e-10. Accuracy: rel.KKT / dx as defined in conic_test_utils.hpp.
//
// Usage: conic_replay [--ipm-tol 1e-8] [--csv out.csv] step_*.bin ...

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "conic_test_utils.hpp"

using namespace conic_test;
namespace scs = cardillo::solver::scs;

namespace {

struct Sample {
    double ms{0}, iters{0}, kkt{0}, dx{0}, setup_ms{0}, scale_updates{0};
    bool ok{false};
};

struct Column {
    std::vector<Sample> s;
    void add(const Sample& x) { s.push_back(x); }
};

double median(std::vector<double> v) {
    if (v.empty()) return std::nan("");
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

double mean(const std::vector<double>& v) {
    double t = 0;
    for (double x : v) t += x;
    return v.empty() ? std::nan("") : t / v.size();
}

}  // namespace

int main(int argc, char** argv) {
    double ipm_tol = 1e-8;
    std::string csv_path;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--ipm-tol") && i + 1 < argc)
            ipm_tol = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc)
            csv_path = argv[++i];
        else
            files.emplace_back(argv[i]);
    }
    if (files.empty()) {
        std::cerr << "usage: conic_replay [--ipm-tol 1e-8] [--csv out.csv] step_*.bin ...\n";
        return 2;
    }

    std::ofstream csv;
    if (!csv_path.empty()) {
        csv.open(csv_path);
        csv << "file,n,m,n_soc,solver,ok,ms,iters,rel_kkt,dx,setup_ms,scale_updates\n";
    }

    std::vector<std::string> order;
    std::map<std::string, Column> cols;
    auto record = [&](const std::string& file, const scs::ProblemDump& d, const std::string& name, const Sample& x) {
        if (!cols.count(name)) order.push_back(name);
        cols[name].add(x);
        if (csv)
            csv << file << ',' << d.P.cols() << ',' << d.A.rows() << ',' << d.dims.soc.size() << ',' << name << ',' << x.ok << ',' << x.ms << ',' << x.iters << ',' << x.kkt << ','
                << x.dx << ',' << x.setup_ms << ',' << x.scale_updates << '\n';
    };

    long n_sum = 0, m_sum = 0, soc_sum = 0;
    for (const auto& file : files) {
        scs::ProblemDump d;
        if (!scs::readDump(file, d)) {
            std::cerr << "cannot read " << file << "\n";
            return 1;
        }
        Problem pr;
        pr.P = d.P;
        pr.q = d.q;
        pr.A = d.A;
        pr.b = d.b;
        pr.dims.zero = d.dims.zero;
        pr.dims.nonneg = d.dims.nonneg;
        pr.dims.soc = d.dims.soc;
        n_sum += pr.P.cols();
        m_sum += pr.A.rows();
        soc_sum += static_cast<long>(pr.dims.soc.size());

        const Result ref = solveClarabel(pr, 1e-10);
        auto ipm = [&](const std::string& name, const Result& r) {
            Sample x;
            x.ok = r.ok;
            x.ms = 1e3 * r.seconds;
            x.iters = std::nan("");
            x.kkt = relKkt(pr, r);
            x.dx = relDx(r, ref);
            record(file, d, name, x);
        };
        ipm("Clarabel", solveClarabel(pr, ipm_tol));
        ipm("QOCO", solveQoco(pr, ipm_tol));
        ipm("ConicXX", solveConicxx(pr, ipm_tol));
#ifdef CARDILLO_HAVE_MOSEK
        {
            cardillo::solver::mosek::Options o;
            o.tol_rel_gap = o.tol_pfeas = o.tol_dfeas = ipm_tol;
            cardillo::solver::mosek::ConicSolver ms(o);
            ipm("MOSEK", solveMosek(ms, pr));
        }
#endif
        for (double eps : {1e-4, 1e-6}) {
            for (bool warm : {false, true}) {
                if (warm && !d.warm) continue;
                scs::Settings st;
                st.eps_abs = st.eps_rel = eps;
                Result w;
                if (warm) {
                    st.scale = d.scale;
                    w.ok = true;
                    w.x = d.x0;
                    w.z = d.y0;
                }
                scs::Info info;
                const Result r = solveScs(pr, st, warm ? &w : nullptr, &info);
                Sample x;
                x.ok = r.ok;
                x.ms = 1e3 * r.seconds;
                x.iters = info.iterations;
                x.kkt = relKkt(pr, r);
                x.dx = relDx(r, ref);
                x.setup_ms = info.setup_ms;
                x.scale_updates = info.scale_updates;
                char name[64];
                std::snprintf(name, sizeof name, "SCS %s %.0e", warm ? "warm" : "cold", eps);
                record(file, d, name, x);
            }
        }
    }

    const double nf = static_cast<double>(files.size());
    std::printf("%zu problems, mean n = %.0f, m = %.0f, SOC blocks = %.0f; IPM tol %.0e; dx reference: Clarabel 1e-10\n\n", files.size(), n_sum / nf, m_sum / nf, soc_sum / nf, ipm_tol);
    std::printf("%-16s %6s %11s %11s %9s %11s %11s %11s %9s\n", "solver", "solved", "median ms", "mean ms", "mean it", "med relKKT", "max relKKT", "med dx", "scale upd");
    for (const auto& name : order) {
        const auto& c = cols[name];
        std::vector<double> ms, it, kkt, dx, su;
        int ok = 0;
        for (const auto& x : c.s) {
            ms.push_back(x.ms);
            if (!std::isnan(x.iters)) it.push_back(x.iters);
            if (x.ok) {
                ++ok;
                kkt.push_back(x.kkt);
                dx.push_back(x.dx);
            }
            su.push_back(x.scale_updates);
        }
        const double max_kkt = kkt.empty() ? std::nan("") : *std::max_element(kkt.begin(), kkt.end());
        std::printf("%-16s %3d/%-3zu %11.3f %11.3f %9.1f %11.1e %11.1e %11.1e %9.1f\n", name.c_str(), ok, c.s.size(), median(ms), mean(ms), it.empty() ? std::nan("") : mean(it), median(kkt), max_kkt,
                    median(dx), mean(su));
    }
    return 0;
}
