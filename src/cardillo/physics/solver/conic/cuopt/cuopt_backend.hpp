#pragma once

#include <memory>
#include <string>

#include "../../../../misc/types.hpp"
#include "cuopt_model.hpp"

namespace cardillo::solver::cuopt {

enum class Presolve { Default, Off, Papilo, Pslp };

/// Settings of the cuOpt barrier method; the defaults of tol_* match ip.tol_abs/ip.tol_rel.
struct Settings {
    real_t tol_abs{1e-8};     // absolute primal, dual and gap tolerances
    real_t tol_rel{1e-8};     // relative primal, dual and gap tolerances
    real_t time_limit{-1};    // per solve [s]; <= 0: no limit
    int iteration_limit{-1};  // < 0: cuOpt's default
    Presolve presolve{Presolve::Default};
    bool deterministic{false};  // deterministic cuDSS factorization (CUOPT_CUDSS_DETERMINISTIC)
    bool verbose{false};        // cuOpt's log on stdout
};

/// Optimal; Inaccurate: an iterate at an iteration, time or numerical limit, or primal feasible only;
/// Failed: no usable solution.
enum class Status { Optimal, Inaccurate, Failed };

const char* toString(Status status);

struct Info {
    Status status{Status::Failed};
    int termination{-1};   // cuOpt termination status, -1 if cuOptSolve returned no solution
    int iterations{-1};    // barrier iterations, -1 if cuOpt does not report them
    real_t solve_time{0};  // cuOpt's own solve time [s]
    real_t objective{0};   // optimal value of the Wolfe dual = -(optimal value of the cone program)
    std::string message;   // reason of a failure
};

/**
 * @brief Solves cone programs (see cuopt_model.hpp) with NVIDIA cuOpt's GPU barrier method.
 *
 * The solver settings are created once; every setup() builds a new cuOpt problem (the C API has no
 * in-place updates) and every solve() starts cold (cuOpt warm-starts only PDLP). Not thread-safe: one
 * instance per thread. The constructor checks the cuOpt version (>= 26.10.1) and its index and float
 * sizes and throws std::runtime_error on a mismatch. Without a usable NVIDIA GPU and driver, solve()
 * returns Status::Failed.
 */
class Solver {
   public:
    explicit Solver(const Settings& settings);
    ~Solver();
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    /// Builds the cuOpt problem of the cone program; only the upper triangle of P is read.
    void setup(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims);
    Status solve();
    /// x and z of the last solve (not Failed), and s = b - Ax, with A and b of the last setup.
    void recover(const CscMatrix& A, const VectorXr& b);

    const VectorXr& x() const;
    const VectorXr& z() const;
    const VectorXr& s() const;
    const Info& info() const;

    /// Writes the problem of the last setup as an MPS file; false on failure.
    bool writeProblem(const std::string& path) const;

    /// Version of the loaded cuOpt library, e.g. "26.10.1".
    static std::string version();

   private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

/**
 * @brief True if cuOpt can solve on this machine: solves a small cone program, which also absorbs
 * the start-up cost of CUDA and cuDSS. Otherwise `reason` (if given) receives the cause.
 */
bool deviceAvailable(std::string* reason = nullptr);

}  // namespace cardillo::solver::cuopt
