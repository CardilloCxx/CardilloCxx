#include <iostream>
#include <memory>
#include <vector>
#include <cstring>

#include "bench/common.hpp"
#include "bench/eigen_dense_ops.hpp"
#include "bench/eigen_sparse_ops.hpp"
#include "bench/block_array_ops.hpp"
#ifdef BENCH_HAVE_BLAS
#include "bench/blas_ops.hpp"
#endif
#ifdef BENCH_HAVE_CUBLAS
#include "bench/cublas_ops.hpp"
#endif
#ifdef BENCH_HAVE_CUSPARSE
#include "bench/cusparse_ops.hpp"
#endif

using bench::IBlockOps;
using bench::Timer;
using bench::Real;

static std::unique_ptr<IBlockOps> make_impl(const std::string& which)
{
    if (which == "dense") return std::make_unique<bench::EigenDenseOps>();
    if (which == "sparse") return std::make_unique<bench::EigenSparseOps>();
    if (which == "block") return std::make_unique<bench::BlockArrayOps>();
#ifdef BENCH_HAVE_BLAS
    if (which == "blas") return std::make_unique<bench::BlasOps>();
#endif
#ifdef BENCH_HAVE_CUBLAS
    if (which == "cublas") return std::make_unique<bench::CuBLASOps>();
#endif
#ifdef BENCH_HAVE_CUSPARSE
    if (which == "cusparse") return std::make_unique<bench::CuSPARSEOps>();
#endif
    std::cerr << "Unknown impl '" << which << "', using block" << std::endl;
    return std::make_unique<bench::BlockArrayOps>();
}

struct Args { std::string impl = "block"; int nb=2000; int nc=4000; int iters=20; unsigned seed=42; double rigidFrac=0.5; bool verifyAll=true; };

static Args parse_args(int argc, char** argv)
{
    Args a;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--impl") && i+1 < argc) a.impl = argv[++i];
        else if (!std::strcmp(argv[i], "--nb") && i+1 < argc) a.nb = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--nc") && i+1 < argc) a.nc = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--iters") && i+1 < argc) a.iters = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--seed") && i+1 < argc) a.seed = (unsigned)std::strtoul(argv[++i], nullptr, 10);
    else if (!std::strcmp(argv[i], "--rigid-frac") && i+1 < argc) a.rigidFrac = std::atof(argv[++i]);
    else if (!std::strcmp(argv[i], "--no-verify")) a.verifyAll = false;
    }
    return a;
}

int main(int argc, char** argv)
{
    auto args = parse_args(argc, argv);
    auto impl = make_impl(args.impl);

    std::cout << "Impl: " << impl->name() << " nb=" << args.nb << " nc=" << args.nc << "\n";

    // Generate mixed DOF per body (3 or 6) and contact pairs
    auto dof = bench::generateDofPerBody(args.nb, args.rigidFrac, args.seed);
    auto contacts = bench::generateContacts(args.nb, args.nc, args.seed);
    impl->buildSystem(dof, contacts, args.seed);

    // Allocate vectors
    int ndof = 0; auto offs = bench::computeOffsets(dof, ndof);
    std::vector<Real> x(ndof), y(3*args.nc), w(3*args.nc), g(ndof), z(ndof);
    bench::fillRandom(x, args.seed+1);
    bench::fillRandom(w, args.seed+2);

    // Warm-up
    impl->mul(x, y);

    // Optional correctness check across all implementations for same seed
    if (args.verifyAll) {
        std::vector<std::unique_ptr<IBlockOps>> impls;
    impls.emplace_back(std::make_unique<bench::BlockArrayOps>());
    impls.emplace_back(std::make_unique<bench::EigenSparseOps>());
    impls.emplace_back(std::make_unique<bench::EigenDenseOps>());
#ifdef BENCH_HAVE_BLAS
    impls.emplace_back(std::make_unique<bench::BlasOps>());
#endif
#ifdef BENCH_HAVE_CUBLAS
    impls.emplace_back(std::make_unique<bench::CuBLASOps>());
#endif
#ifdef BENCH_HAVE_CUSPARSE
    impls.emplace_back(std::make_unique<bench::CuSPARSEOps>());
#endif
        std::vector<std::vector<Real>> resultsY, resultsZ;
        for (auto& p : impls) {
            p->buildSystem(dof, contacts, args.seed);
            std::vector<Real> y2(3*args.nc), z2(ndof);
            p->mul(x, y2);
            p->normalMatrixMul(x, z2);
            resultsY.push_back(std::move(y2));
            resultsZ.push_back(std::move(z2));
        }
        auto almostEqual = [](const std::vector<Real>& a, const std::vector<Real>& b){
            if (a.size()!=b.size()) return false;
            double maxAbs=0.0, maxRel=0.0; bool ok=true; const double eps=1e-9;
            for (size_t i=0;i<a.size();++i){ double da=std::abs((double)a[i]-b[i]); double m=std::max(1.0,std::abs((double)a[i])+std::abs((double)b[i])); maxAbs=std::max(maxAbs,da); maxRel=std::max(maxRel,da/m); if(da>1e-6 && da/m>1e-6) ok=false; }
            if(!ok) std::cerr<<"Verification failed: maxAbs="<<maxAbs<<" maxRel="<<maxRel<<"\n";
            return ok;
        };
        bool y_ok = almostEqual(resultsY[0], resultsY[1]) && almostEqual(resultsY[0], resultsY[2]);
        bool z_ok = almostEqual(resultsZ[0], resultsZ[1]) && almostEqual(resultsZ[0], resultsZ[2]);
        std::cout << "Verify A*x equal: " << (y_ok?"OK":"FAIL") << ", AtAx equal: " << (z_ok?"OK":"FAIL") << "\n";
    }

    // Benchmark A*x
    Timer t; t.tic();
    for (int it = 0; it < args.iters; ++it) impl->mul(x, y);
    double t_mul = t.tocMs();

    // Benchmark rowDot over all rows
    t.tic();
    Real acc = 0;
    for (int it = 0; it < args.iters; ++it) {
        for (int k = 0; k < args.nc; ++k) acc += impl->rowDot((std::size_t)k, x);
    }
    double t_row = t.tocMs();

    // Benchmark A^T * w accumulation
    std::fill(g.begin(), g.end(), Real(0));
    t.tic();
    for (int it = 0; it < args.iters; ++it) impl->mulTransposeAcc(w, g);
    double t_at = t.tocMs();

    // Benchmark normal matrix multiply
    t.tic();
    for (int it = 0; it < args.iters; ++it) impl->normalMatrixMul(x, z);
    double t_normal = t.tocMs();

    // Basic checksums to avoid dead-code elimination
    auto checksum = [&](const std::vector<Real>& v){ double s=0; for (auto& e: v) s += e; return s; };
    std::cout << "Timing (" << args.iters << " iters):\n";
    std::cout << "  A*x:         " << t_mul << " ms, sum(y)=" << checksum(y) << "\n";
    std::cout << "  rowDot sum:  " << t_row << " ms, acc=" << acc << "\n";
    std::cout << "  At*w (acc):  " << t_at << " ms, sum(g)=" << checksum(g) << "\n";
    std::cout << "  At*A*x:      " << t_normal << " ms, sum(z)=" << checksum(z) << "\n";

    return 0;
}
