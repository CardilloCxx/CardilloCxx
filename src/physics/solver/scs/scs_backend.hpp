#pragma once

#include <string>
#include <vector>

#include <Eigen/SparseCore>

#include "../../../misc/types.hpp"

// Opaque SCS workspace (defined in SCS's scs_work.h); scs.h itself is only included by the .cpp.
struct SCS_WORK;

namespace cardillo::solver::scs {

using SpMat = Eigen::SparseMatrix<real_t, Eigen::ColMajor, int>;

/**
 * @brief Cone composition of the rows of the stacked constraint matrix, in row order: `zero`
 * equality rows, then one `nonneg`-dimensional nonnegative orthant, then one second-order cone per
 * entry of `soc` (scalar component first). This is the row layout ClarabelAssembler produces, and it
 * is also the order SCS requires (zero, positive, ..., second-order), so no row permutation is needed.
 */
struct ConeDims {
    int zero{0};
    int nonneg{0};
    std::vector<int> soc;

    int rows() const;
};

/**
 * @brief SCS parameters. Negative numeric values leave SCS's own default (scs_set_default_settings)
 * in place.
 */
struct Settings {
    double eps_abs{-1};
    double eps_rel{-1};
    double eps_infeas{-1};
    int max_iters{-1};
    double time_limit_secs{-1};
    double alpha{-1};
    double scale{-1};
    double rho_x{-1};
    int acceleration_lookback{-1};
    int acceleration_interval{-1};
    bool adaptive_scale{true};
    bool normalize{true};
    bool verbose{false};
};

enum class Status { Solved, SolvedInaccurate, Infeasible, Unbounded, Failed, Interrupted, Unknown };

const char* toString(Status s);
/// Solved or SolvedInaccurate: the iterate is a usable (approximate) primal-dual solution.
inline bool isAcceptable(Status s) { return s == Status::Solved || s == Status::SolvedInaccurate; }

/// Per-solve statistics as reported by SCS (times in milliseconds).
struct Info {
    Status status{Status::Unknown};
    int iterations{0};
    double setup_ms{0};    // scs_init: validation, normalization, AMD + symbolic + numeric factorization
    double solve_ms{0};    // scs_solve
    double lin_sys_ms{0};  // part of solve_ms spent in the KKT solves (incl. refactorizations on scale updates)
    double cone_ms{0};     // part of solve_ms spent in cone projections
    double accel_ms{0};    // part of solve_ms spent in Anderson acceleration
    int accel_accepted{0};
    int accel_rejected{0};
    int scale_updates{0};
    double scale{0};  // final (adapted) dual scale
    double res_pri{0};
    double res_dual{0};
    double gap{0};
    double pobj{0};
};

/**
 * @brief Thin owning wrapper around SCS's three-phase API (scs_init / scs_solve / scs_finish) with
 * the direct (AMD + QDLDL) linear-system backend, for problems in Clarabel's form
 *
 *   min 1/2 x'Px + q'x   s.t.  A x + s = b,  s in K,
 *
 * with the dual y satisfying P x + q + A'y = 0, y in K*. SCS keeps deep copies of the data, so the
 * inputs only have to live for the duration of setup().
 */
class Solver {
   public:
    Solver();
    ~Solver();
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    /// scs_init on new data (any previous workspace is released). P may be full or upper-triangular;
    /// only its upper triangle is passed on. Returns false if SCS rejects the data.
    bool setup(const SpMat& P, const VectorXr& q, const SpMat& A, const VectorXr& b, const ConeDims& dims, const Settings& settings);

    /// Mutable warm-start buffers (sized by setup()). Filled by the caller before solve(true); after
    /// solve() they hold the solution.
    VectorXr& x() { return m_x; }
    VectorXr& y() { return m_y; }
    VectorXr& s() { return m_s; }
    const VectorXr& x() const { return m_x; }
    const VectorXr& y() const { return m_y; }
    const VectorXr& s() const { return m_s; }

    Status solve(bool warm_start);

    const Info& info() const { return m_info; }
    int n() const { return static_cast<int>(m_x.size()); }
    int m() const { return static_cast<int>(m_y.size()); }

   private:
    void release();

    SCS_WORK* m_work{nullptr};
    VectorXr m_x, m_y, m_s;
    Info m_info;
};

/**
 * @brief One step's cone program plus the warm start the solver built for it, for replaying exactly
 * the same data through other solvers offline (scs.dump_dir, tests/conic_replay.cpp). Binary,
 * native endianness, not meant as an exchange format.
 */
struct ProblemDump {
    SpMat P;  // as assembled (ClarabelAssembler: diagonal, i.e. upper == full)
    VectorXr q;
    SpMat A;
    VectorXr b;
    ConeDims dims;
    bool warm{false};
    double scale{-1};  // initial dual scale used for this solve (-1: SCS default)
    VectorXr x0, y0, s0;
};

bool writeDump(const std::string& path, const ProblemDump& d);
bool readDump(const std::string& path, ProblemDump& d);

/// Euclidean projection of the cone part (rows after `dims.zero`) of `v` onto K; the zero-cone part
/// is left untouched. K is self-dual apart from the zero cone, so this projects onto K* as well.
void projectOntoCones(VectorXr& v, const ConeDims& dims);

}  // namespace cardillo::solver::scs
