#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <mosek.h>

#include "../../misc/types.hpp"

namespace cardillo::solver::mosek {

// ---------------------------------------------------------------------------------------------
// Error handling / RAII
// ---------------------------------------------------------------------------------------------

/**
 * @brief Throws std::runtime_error if `r != MSK_RES_OK`. The message names the failed `operation`,
 * MOSEK's symbolic code name and description and, when a task is given, MOSEK's last error message
 * for that task, e.g. "MOSEK error during MSK_putafefentrylist: MSK_RES_ERR_INDEX (1235): ...".
 */
void check(MSKrescodee r, const char* operation, MSKtask_t task = nullptr);

/// Owning wrapper of an MSKenv_t (one per solver; it also holds the checked-out license).
class Environment {
   public:
    Environment();
    ~Environment();
    Environment(const Environment&) = delete;
    Environment& operator=(const Environment&) = delete;

    MSKenv_t get() const { return m_env; }

   private:
    MSKenv_t m_env{nullptr};
};

/// Owning wrapper of an MSKtask_t. reset() destroys the current task and creates an empty one.
class Task {
   public:
    Task() = default;
    ~Task();
    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;

    void reset(const Environment& env, MSKint32t maxnumcon, MSKint32t maxnumvar);
    MSKtask_t get() const { return m_task; }
    explicit operator bool() const { return m_task != nullptr; }

   private:
    void destroy();
    MSKtask_t m_task{nullptr};
};

// ---------------------------------------------------------------------------------------------
// Problem description, options, results
// ---------------------------------------------------------------------------------------------

/**
 * @brief Cone composition of the rows of the stacked constraint matrix, in row order: `zero`
 * equality rows, then one `nonneg`-dimensional nonnegative orthant, then one second-order cone per
 * entry of `soc`. This is the same row layout as Clarabel's ZeroCone/NonnegativeCone/SecondOrderCone
 * list and QOCO's (A; G) split with `p = zero`, `l = nonneg`, `q = soc`.
 */
struct ConeDims {
    int zero{0};
    int nonneg{0};
    std::vector<int> soc;

    int rows() const;
    bool operator==(const ConeDims& o) const { return zero == o.zero && nonneg == o.nonneg && soc == o.soc; }
};

/**
 * @brief MOSEK parameters. Negative values (or presolve = Default) leave MOSEK's own default in
 * place; nothing is tuned by default.
 */
struct Options {
    enum class Presolve { Default, Off, On, Free };
    /// How 1/2 x'Px enters the MOSEK model (see ConicSolver). Native (MSK_putqobj) is only
    /// possible without cone rows; Separable needs a diagonal P. Auto picks Native without cone
    /// rows, else Separable for a diagonal P and Single otherwise.
    enum class QuadModel { Auto, Native, Single, Separable };

    int num_threads{0};           // MSK_IPAR_NUM_THREADS (0 = let MOSEK decide)
    double tol_rel_gap{-1.0};     // MSK_DPAR_INTPNT_CO_TOL_REL_GAP
    double tol_pfeas{-1.0};       // MSK_DPAR_INTPNT_CO_TOL_PFEAS
    double tol_dfeas{-1.0};       // MSK_DPAR_INTPNT_CO_TOL_DFEAS
    int max_iterations{-1};       // MSK_IPAR_INTPNT_MAX_ITERATIONS
    double max_time{-1.0};        // MSK_DPAR_OPTIMIZER_MAX_TIME [s]
    Presolve presolve{Presolve::Default};  // MSK_IPAR_PRESOLVE_USE
    QuadModel quad_model{QuadModel::Auto};
    bool scaling{true};           // diagonal pre-scaling, see ConicSolver
    bool verbose{false};          // forward MOSEK's log stream to std::cout
    bool equality_duals{true};    // also retrieve the duals of the zero-cone (equality) rows
};

enum class Status {
    Optimal,
    NearOptimal,  // primal and dual feasible, but the gap did not reach its tolerance (e.g. stall)
    PrimalInfeasible,
    DualInfeasible,
    IterationLimit,
    TimeLimit,
    NumericalFailure,
    Unknown,
};

const char* toString(Status s);
inline bool isAcceptable(Status s) { return s == Status::Optimal || s == Status::NearOptimal; }

// ---------------------------------------------------------------------------------------------
// Conic QP solver
// ---------------------------------------------------------------------------------------------

/**
 * @brief Sparse convex conic QP solver on top of the MOSEK Optimizer API (mosek.h, no Fusion).
 *
 * Problem (identical to Clarabel's form, i.e. QOCO's `Ax = b, Gx + s = h, s in K` with A and G
 * stacked row-wise):
 *
 *     min  1/2 x^T P x + q^T x
 *     s.t. A x + s = b,   s in {0}^zero x R_+^nonneg x Q^soc[0] x Q^soc[1] x ...
 *
 * with Q^k = { s in R^k : s_0 >= ||s_{1:}|| } (the same ordering QOCO, Clarabel and MOSEK use, so
 * no cone entries are permuted). Only the lower triangle of P is read; P must be PSD.
 *
 * MOSEK formulation (variables [x; t]):
 *   - zero rows:       ordinary linear constraints, bounds l_c = u_c = b (MSK_BK_FX).
 *   - cone rows:       affine conic constraints  F x + g in K  with  F = -A,  g = b,
 *                      i.e.  b - A x = s in K.  The sign flip is applied while copying values;
 *                      no negated matrix is formed.
 *   - quadratic term:  without cone rows (pure equality-constrained QP) P's lower triangle is
 *                      passed natively with MSK_putqobj (MOSEK's objective also carries the 1/2).
 *                      MOSEK rejects a quadratic objective combined with conic constraints
 *                      (MSK_RES_ERR_MIXED_CONIC_AND_NL), so otherwise 1/2 x^T P x is modelled
 *                      with rotated-cone epigraphs (Q_r^k = {2 u0 u1 >= ||u_{2:}||^2, u0,u1 >= 0}),
 *                      where L^T L = P:
 *                      - Separable (default for diagonal P, the Cardillo case): one epigraph
 *                        variable t_i per row of L = diag(sqrt(P_ii)), (t_i, 1, L_i x) in Q_r^3,
 *                        objective q^T x + sum t_i. Each cone couples only one variable, which
 *                        keeps the interior-point iteration well conditioned for large n.
 *                      - Single: one variable t with (t, 1, L x) in Q_r^{r+2}, objective
 *                        q^T x + t. Used for non-diagonal P, where L comes from a sparse LDL^T
 *                        factorization of P (P must then be positive definite).
 *
 * Scaling (Options::scaling, on by default): MOSEK's own scaling does not cope with the mass/
 * compliance ranges of Cardillo problems (P_ii spanning ~1e-10..1e7 stalls it), so the problem is
 * passed in the variables x = D xh, D_jj = 1/sqrt(P_jj) (1 where P_jj = 0), with every zero-cone row
 * and every cone block (one common positive factor per block, which preserves the cone) scaled by
 * the inverse of its largest |(A D)_ij|. This is applied while copying the coefficients, and x and
 * the duals are mapped back in extract() (x = D xh, z = E zh), so callers never see it.
 *
 * Dual sign convention of the returned duals (matches Clarabel): P x + q + A^T z = 0, with the cone
 * part of z in K* (= K, all cones here are self-dual).
 *
 * The task is persistent: setup() builds the structure (variables, constraints, AFEs, domains,
 * ACCs); update() only overwrites numerical values and is valid as long as structureMatches()
 * returns true for the new data (same dimensions, cone dims, and sparsity patterns of P and A).
 * All index/value buffers are kept between calls, so the update path does not allocate.
 *
 * Warm starts: MOSEK's interior-point optimizer does not accept an initial point, so none is
 * provided here.
 */
class ConicSolver {
   public:
    explicit ConicSolver(const Options& options = {});

    void setOptions(const Options& options);
    const Options& options() const { return m_opts; }

    bool hasStructure() const { return static_cast<bool>(m_task); }
    bool structureMatches(const SparseMatrix<Eigen::ColMajor>& P, const SparseMatrix<Eigen::ColMajor>& A, const ConeDims& dims) const;

    /// Builds a fresh task for the given problem (also loads all numerical values).
    void setup(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims);
    /// Overwrites the numerical values only. Requires structureMatches(P, A, dims).
    void update(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b);
    /// setup() or update(), whichever applies. Returns true if the task was (re)built.
    bool load(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims);

    Status optimize();
    /// Copies x and z (and the equality duals if enabled) out of the task.
    void extract();

    /// Convenience: load() + optimize() + extract() (extract only if a solution is available).
    Status solve(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims);

    // Results of the last optimize()/extract().
    Status status() const { return m_status; }
    Eigen::VectorBlock<const VectorXr> x() const { return m_xx.head(m_n); }
    /// Duals of the cone (non-zero-cone) rows, in row order.
    Eigen::VectorBlock<const VectorXr> coneDuals() const { return m_doty.head(m_m); }
    /// Duals of the zero-cone rows (empty unless Options::equality_duals).
    const VectorXr& equalityDuals() const { return m_y; }
    double objective() const { return m_objective; }
    int iterations() const { return m_iterations; }
    double optimizerTime() const { return m_optimizer_time; }
    MSKrescodee terminationCode() const { return m_trm; }

   private:
    void applyOptions();
    void computeQuadFactor(const SparseMatrix<Eigen::ColMajor>& P, bool structural);
    void computeScaling(const SparseMatrix<Eigen::ColMajor>& P, const SparseMatrix<Eigen::ColMajor>& A);
    void fillValues(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b);
    void putValues();

    Options m_opts;
    Environment m_env;
    Task m_task;

    // Structure of the loaded problem.
    int m_n{0};        // number of x variables
    int m_p{0};        // equality rows
    int m_m{0};        // cone rows
    int m_r{0};        // rows of the quadratic factor L
    ConeDims m_dims;
    std::vector<int> m_A_outer, m_A_inner, m_P_outer, m_P_inner;
    bool m_P_diagonal{true};
    bool m_native{false};     // quadratic term via MSK_putqobj, see Options::QuadModel
    bool m_separable{false};  // epigraph layout, see Options::QuadModel
    Options::QuadModel m_built_quad_model{Options::QuadModel::Auto};  // option the task was built with

    // Native quadratic objective: lower-triangle entries of P.
    std::vector<MSKint32t> m_q_subi, m_q_subj;
    std::vector<MSKint64t> m_q_src;  // index into P.valuePtr()
    std::vector<MSKrealt> m_q_val;

    // Cached index arrays (structure) and value arrays (overwritten in place on update).
    std::vector<MSKint32t> m_a_subi, m_a_subj;  // linear (equality) constraint entries
    std::vector<MSKrealt> m_a_val;
    std::vector<MSKint64t> m_f_afe;  // AFE entries: cone rows of -A, then (t) and L
    std::vector<MSKint32t> m_f_var;
    std::vector<MSKrealt> m_f_val;
    std::vector<MSKint64t> m_lin_src;   // index into A.valuePtr() of each linear entry
    std::vector<MSKint64t> m_cone_src;  // index into A.valuePtr() of each F entry taken from A
    std::vector<MSKboundkeye> m_bkc;
    std::vector<MSKrealt> m_blc;
    std::vector<MSKrealt> m_g;  // AFE constant terms
    std::vector<MSKrealt> m_c;  // (scaled) linear objective of x
    VectorXr m_D;               // variable scaling, x = D xh
    VectorXr m_E;               // row scaling of A's rows

    // Quadratic factor L (rows m_r, cols m_n) as triplets; pattern fixed for a fixed P pattern.
    // Non-diagonal P only: LDL^T factorization (symbolic analysis cached with the structure).
    Eigen::SimplicialLDLT<SparseMatrix<Eigen::ColMajor>, Eigen::Lower> m_ldlt;
    std::vector<MSKint32t> m_l_row, m_l_col;
    std::vector<MSKrealt> m_l_val;

    // Results
    Status m_status{Status::Unknown};
    MSKrescodee m_trm{MSK_RES_OK};
    VectorXr m_xx, m_doty, m_y;
    double m_objective{0.0};
    int m_iterations{0};
    double m_optimizer_time{0.0};
};

}  // namespace cardillo::solver::mosek
