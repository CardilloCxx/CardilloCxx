#include "mosek_backend.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace cardillo::solver::mosek {

namespace {

void MSKAPI logToStdout(MSKuserhandle_t, const char* str) { std::cout << str; }

template <typename Index>
bool samePattern(const SparseMatrix<Eigen::ColMajor>& M, const std::vector<Index>& outer, const std::vector<Index>& inner) {
    const auto nouter = static_cast<std::size_t>(M.outerSize()) + 1;
    const auto nnz = static_cast<std::size_t>(M.nonZeros());
    return outer.size() == nouter && inner.size() == nnz && std::equal(outer.begin(), outer.end(), M.outerIndexPtr()) &&
           std::equal(inner.begin(), inner.end(), M.innerIndexPtr());
}

void storePattern(const SparseMatrix<Eigen::ColMajor>& M, std::vector<int>& outer, std::vector<int>& inner) {
    outer.assign(M.outerIndexPtr(), M.outerIndexPtr() + M.outerSize() + 1);
    inner.assign(M.innerIndexPtr(), M.innerIndexPtr() + M.nonZeros());
}

void requireCompressed(const SparseMatrix<Eigen::ColMajor>& M, const char* name) {
    if (!M.isCompressed()) throw std::invalid_argument(std::string("mosek::ConicSolver: matrix ") + name + " must be in compressed storage");
}

}  // namespace

// ---------------------------------------------------------------------------------------------

void check(MSKrescodee r, const char* operation, MSKtask_t task) {
    if (r == MSK_RES_OK) return;

    std::string msg = std::string("MOSEK error during ") + operation + ": ";
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    if (MSK_getcodedesc(r, symname, desc) == MSK_RES_OK)
        msg += std::string(symname) + " (" + std::to_string(static_cast<int>(r)) + "): " + desc;
    else
        msg += "code " + std::to_string(static_cast<int>(r));

    if (task) {
        MSKrescodee last = MSK_RES_OK;
        MSKint32t len = 0;
        char lastmsg[MSK_MAX_STR_LEN];
        if (MSK_getlasterror(task, &last, MSK_MAX_STR_LEN, &len, lastmsg) == MSK_RES_OK && len > 0) msg += std::string("\n  ") + lastmsg;
    }
    throw std::runtime_error(msg);
}

Environment::Environment() { check(MSK_makeenv(&m_env, nullptr), "MSK_makeenv"); }

Environment::~Environment() {
    if (m_env) MSK_deleteenv(&m_env);
}

Task::~Task() { destroy(); }

void Task::destroy() {
    if (m_task) MSK_deletetask(&m_task);
    m_task = nullptr;
}

void Task::reset(const Environment& env, MSKint32t maxnumcon, MSKint32t maxnumvar) {
    destroy();
    check(MSK_maketask(env.get(), maxnumcon, maxnumvar, &m_task), "MSK_maketask");
}

int ConeDims::rows() const {
    int r = zero + nonneg;
    for (int d : soc) r += d;
    return r;
}

const char* toString(Status s) {
    switch (s) {
        case Status::Optimal:
            return "Optimal";
        case Status::NearOptimal:
            return "Near Optimal";
        case Status::PrimalInfeasible:
            return "Primal Infeasible";
        case Status::DualInfeasible:
            return "Dual Infeasible";
        case Status::IterationLimit:
            return "Iteration Limit";
        case Status::TimeLimit:
            return "Time Limit";
        case Status::NumericalFailure:
            return "Numerical Failure";
        case Status::Unknown:
        default:
            return "Unknown";
    }
}

// ---------------------------------------------------------------------------------------------

ConicSolver::ConicSolver(const Options& options) : m_opts(options) {}

void ConicSolver::setOptions(const Options& options) {
    m_opts = options;
    if (m_task) applyOptions();
}

void ConicSolver::applyOptions() {
    MSKtask_t task = m_task.get();
    check(MSK_putintparam(task, MSK_IPAR_NUM_THREADS, std::max(0, m_opts.num_threads)), "MSK_putintparam(NUM_THREADS)", task);
    if (m_opts.tol_rel_gap > 0) check(MSK_putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, m_opts.tol_rel_gap), "MSK_putdouparam(INTPNT_CO_TOL_REL_GAP)", task);
    if (m_opts.tol_pfeas > 0) check(MSK_putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_PFEAS, m_opts.tol_pfeas), "MSK_putdouparam(INTPNT_CO_TOL_PFEAS)", task);
    if (m_opts.tol_dfeas > 0) check(MSK_putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_DFEAS, m_opts.tol_dfeas), "MSK_putdouparam(INTPNT_CO_TOL_DFEAS)", task);
    if (m_opts.max_iterations > 0) check(MSK_putintparam(task, MSK_IPAR_INTPNT_MAX_ITERATIONS, m_opts.max_iterations), "MSK_putintparam(INTPNT_MAX_ITERATIONS)", task);
    if (m_opts.max_time > 0) check(MSK_putdouparam(task, MSK_DPAR_OPTIMIZER_MAX_TIME, m_opts.max_time), "MSK_putdouparam(OPTIMIZER_MAX_TIME)", task);
    if (m_opts.presolve != Options::Presolve::Default) {
        const MSKpresolvemodee mode = (m_opts.presolve == Options::Presolve::Off)  ? MSK_PRESOLVE_MODE_OFF
                                      : (m_opts.presolve == Options::Presolve::On) ? MSK_PRESOLVE_MODE_ON
                                                                                   : MSK_PRESOLVE_MODE_FREE;
        check(MSK_putintparam(task, MSK_IPAR_PRESOLVE_USE, mode), "MSK_putintparam(PRESOLVE_USE)", task);
    }
    // The problem always has ACCs unless it is a pure LP; pin the interior-point optimizer so the
    // interior-point solution (MSK_SOL_ITR) is the one that is always defined, without crossover.
    const bool conic = m_m > 0 || m_r > 0;
    check(MSK_putintparam(task, MSK_IPAR_OPTIMIZER, conic ? MSK_OPTIMIZER_CONIC : MSK_OPTIMIZER_INTPNT), "MSK_putintparam(OPTIMIZER)", task);
    check(MSK_putintparam(task, MSK_IPAR_INTPNT_BASIS, MSK_BI_NEVER), "MSK_putintparam(INTPNT_BASIS)", task);
}

bool ConicSolver::structureMatches(const SparseMatrix<Eigen::ColMajor>& P, const SparseMatrix<Eigen::ColMajor>& A, const ConeDims& dims) const {
    if (!m_task) return false;
    if (P.rows() != m_n || P.cols() != m_n || A.cols() != m_n || A.rows() != m_p + m_m) return false;
    if (!(dims == m_dims) || m_built_quad_model != m_opts.quad_model) return false;
    if (!P.isCompressed() || !A.isCompressed()) return false;
    return samePattern(A, m_A_outer, m_A_inner) && samePattern(P, m_P_outer, m_P_inner);
}

bool ConicSolver::load(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims) {
    if (structureMatches(P, A, dims)) {
        update(P, q, A, b);
        return false;
    }
    setup(P, q, A, b, dims);
    return true;
}

Status ConicSolver::solve(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims) {
    load(P, q, A, b, dims);
    const Status s = optimize();
    if (s != Status::PrimalInfeasible && s != Status::DualInfeasible) extract();
    return s;
}

// Builds the triplets of the quadratic factor L with L^T L = P (only P's lower triangle is read).
// `structural`: (re)compute the index arrays; otherwise only the values are overwritten.
void ConicSolver::computeQuadFactor(const SparseMatrix<Eigen::ColMajor>& P, bool structural) {
    const int* outer = P.outerIndexPtr();
    const int* inner = P.innerIndexPtr();
    const double* val = P.valuePtr();

    if (structural) {
        m_P_diagonal = true;
        for (int j = 0; j < m_n && m_P_diagonal; ++j)
            for (int k = outer[j]; k < outer[j + 1]; ++k)
                if (inner[k] != j) {
                    m_P_diagonal = false;
                    break;
                }
    }

#ifndef NDEBUG
    if (!m_P_diagonal) {
        const double asym = (SparseMatrix<Eigen::ColMajor>(P.transpose()) - P).norm();
        if (asym > 1e-10 * std::max(1.0, P.norm())) std::cerr << "Warning: mosek::ConicSolver: P is not symmetric (||P^T - P|| = " << asym << "); only its lower triangle is used.\n";
    }
#endif

    if (m_P_diagonal) {
        // L = diag(sqrt(P_jj)) restricted to the stored diagonal entries.
        if (structural) {
            m_l_row.clear();
            m_l_col.clear();
            for (int j = 0; j < m_n; ++j)
                if (outer[j + 1] > outer[j]) {
                    m_l_row.push_back(static_cast<MSKint32t>(m_l_row.size()));
                    m_l_col.push_back(j);
                }
            m_l_val.resize(m_l_row.size());
        }
        for (std::size_t e = 0; e < m_l_col.size(); ++e) {
            const double d = val[outer[m_l_col[e]]];
            if (d < 0) throw std::invalid_argument("mosek::ConicSolver: P has a negative diagonal entry (objective is not convex)");
            m_l_val[e] = std::sqrt(d);
        }
        m_r = static_cast<int>(m_l_row.size());
        return;
    }

    // General P = Pi^T L D L^T Pi (Eigen factorizes Pi P Pi^T with y = Pi x, y_{sigma(j)} = x_j),
    // so x^T P x = || sqrt(D) L^T y ||^2. Row i of sqrt(D) L^T has sqrt(D_i) at y_i and
    // sqrt(D_i) L_ki at y_k (k > i); y_k is x_{sigma^-1(k)}.
    if (structural) m_ldlt.analyzePattern(P);
    m_ldlt.factorize(P);
    if (m_ldlt.info() != Eigen::Success)
        throw std::invalid_argument("mosek::ConicSolver: LDL^T factorization of non-diagonal P failed (P must be positive definite)");

    const auto& L = m_ldlt.matrixL().nestedExpression();
    const auto& D = m_ldlt.vectorD();
    const auto& pinv = m_ldlt.permutationPinv().indices();

    std::size_t e = 0;
    auto emit = [&](int row, int ycol, double v) {
        if (structural) {
            m_l_row.push_back(row);
            m_l_col.push_back(pinv[ycol]);
            m_l_val.push_back(v);
        } else {
            if (e >= m_l_row.size() || m_l_row[e] != row || m_l_col[e] != pinv[ycol])
                throw std::logic_error("mosek::ConicSolver: LDL^T pattern changed although P's pattern did not");
            m_l_val[e] = v;
        }
        ++e;
    };
    if (structural) {
        m_l_row.clear();
        m_l_col.clear();
        m_l_val.clear();
    }
    for (int i = 0; i < m_n; ++i) {
        if (D[i] < 0) throw std::invalid_argument("mosek::ConicSolver: P is not positive semidefinite (negative LDL^T pivot)");
        const double sd = std::sqrt(D[i]);
        emit(i, i, sd);
        for (SparseMatrix<Eigen::ColMajor>::InnerIterator it(L, i); it; ++it)
            if (it.row() > i) emit(i, static_cast<int>(it.row()), sd * it.value());
    }
    m_r = m_n;
}

void ConicSolver::setup(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b, const ConeDims& dims) {
    const int n = static_cast<int>(P.cols());
    if (P.rows() != n || A.cols() != n || q.size() != n) throw std::invalid_argument("mosek::ConicSolver: inconsistent dimensions of P, q and A");
    if (A.rows() != dims.rows() || b.size() != A.rows()) throw std::invalid_argument("mosek::ConicSolver: rows of A/b do not match the cone dimensions");
    if (dims.zero < 0 || dims.nonneg < 0 || std::any_of(dims.soc.begin(), dims.soc.end(), [](int d) { return d < 1; }))
        throw std::invalid_argument("mosek::ConicSolver: invalid cone dimensions");
    requireCompressed(P, "P");
    requireCompressed(A, "A");

    // Invalidate first so a throw below leaves the solver without a (half-built) structure.
    m_task.reset(m_env, dims.zero, n + 1);

    m_n = n;
    m_p = dims.zero;
    m_m = A.rows() - dims.zero;
    m_dims = dims;
    storePattern(A, m_A_outer, m_A_inner);
    storePattern(P, m_P_outer, m_P_inner);

    const Options::QuadModel model = m_opts.quad_model;
    m_built_quad_model = model;
    m_native = model == Options::QuadModel::Native || (model == Options::QuadModel::Auto && m_m == 0);
    if (m_native && m_m > 0) throw std::invalid_argument("mosek::ConicSolver: a native quadratic objective cannot be combined with cone rows");
    if (m_native) {
        // Lower triangle of P as triplets; values are refilled by fillQuadValues().
        m_q_subi.clear();
        m_q_subj.clear();
        m_q_src.clear();
        const int* pouter = P.outerIndexPtr();
        const int* pinner = P.innerIndexPtr();
        for (int j = 0; j < m_n; ++j)
            for (int k = pouter[j]; k < pouter[j + 1]; ++k)
                if (pinner[k] >= j) {
                    m_q_subi.push_back(pinner[k]);
                    m_q_subj.push_back(j);
                    m_q_src.push_back(k);
                }
        m_q_val.resize(m_q_src.size());
        m_r = 0;
        m_l_row.clear();
        m_l_col.clear();
        m_l_val.clear();
    } else {
        m_q_subi.clear();
        m_q_subj.clear();
        m_q_src.clear();
        m_q_val.clear();
        computeQuadFactor(P, true);
    }
    const bool has_quad = m_r > 0;
    if (!m_native && model == Options::QuadModel::Separable && !m_P_diagonal)
        throw std::invalid_argument("mosek::ConicSolver: the separable quadratic epigraph requires a diagonal P");
    m_separable = m_P_diagonal && model != Options::QuadModel::Single;
    const int nquad = has_quad ? (m_separable ? m_r : 1) : 0;  // epigraph variables t
    const int nvar = m_n + nquad;                                // x, then t
    const MSKint64t numafe = static_cast<MSKint64t>(m_m) + (has_quad ? (m_separable ? 3 * m_r : m_r + 2) : 0);

    // Split A's entries by row: zero-cone rows -> linear constraints, cone rows -> F = -A.
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    m_a_subi.clear();
    m_a_subj.clear();
    m_lin_src.clear();
    m_f_afe.clear();
    m_f_var.clear();
    m_cone_src.clear();
    for (int j = 0; j < m_n; ++j) {
        for (int k = outer[j]; k < outer[j + 1]; ++k) {
            const int i = inner[k];
            if (i < m_p) {
                m_a_subi.push_back(i);
                m_a_subj.push_back(j);
                m_lin_src.push_back(k);
            } else {
                m_f_afe.push_back(i - m_p);
                m_f_var.push_back(j);
                m_cone_src.push_back(k);
            }
        }
    }
    if (has_quad && m_separable) {
        // AFE rows m+3e, m+3e+1, m+3e+2 hold (t_e, 1, L_e x) of rotated cone e (L diagonal).
        for (std::size_t e = 0; e < m_l_row.size(); ++e) {
            const MSKint64t base = static_cast<MSKint64t>(m_m) + 3 * static_cast<MSKint64t>(e);
            m_f_afe.push_back(base);
            m_f_var.push_back(m_n + static_cast<int>(e));
            m_f_afe.push_back(base + 2);
            m_f_var.push_back(m_l_col[e]);
        }
    } else if (has_quad) {
        // AFE rows m, m+1, m+2.. hold (t, 1, L x) of the rotated cone.
        m_f_afe.push_back(m_m);
        m_f_var.push_back(m_n);
        for (std::size_t e = 0; e < m_l_row.size(); ++e) {
            m_f_afe.push_back(static_cast<MSKint64t>(m_m) + 2 + m_l_row[e]);
            m_f_var.push_back(m_l_col[e]);
        }
    }
    m_a_val.resize(m_a_subi.size());
    m_f_val.resize(m_f_afe.size());
    m_bkc.assign(static_cast<std::size_t>(m_p), MSK_BK_FX);
    m_blc.resize(static_cast<std::size_t>(m_p));
    m_g.assign(static_cast<std::size_t>(numafe), 0.0);
    if (has_quad && m_separable)
        for (int e = 0; e < m_r; ++e) m_g[static_cast<std::size_t>(m_m) + 3 * static_cast<std::size_t>(e) + 1] = 1.0;
    else if (has_quad)
        m_g[static_cast<std::size_t>(m_m) + 1] = 1.0;

    MSKtask_t task = m_task.get();
    if (m_opts.verbose) check(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, nullptr, logToStdout), "MSK_linkfunctotaskstream", task);
    applyOptions();

    check(MSK_appendvars(task, nvar), "MSK_appendvars", task);
    check(MSK_putvarboundsliceconst(task, 0, nvar, MSK_BK_FR, -MSK_INFINITY, MSK_INFINITY), "MSK_putvarboundsliceconst", task);
    check(MSK_appendcons(task, m_p), "MSK_appendcons", task);
    check(MSK_appendafes(task, numafe), "MSK_appendafes", task);
    if (has_quad) check(MSK_putcslice(task, m_n, nvar, std::vector<MSKrealt>(static_cast<std::size_t>(nquad), 1.0).data()), "MSK_putcslice(t)", task);

    // Domains (one per distinct cone dimension, shared by all ACCs of that dimension) and ACCs,
    // laid out contiguously over AFEs 0..numafe-1 in row order.
    std::vector<MSKint64t> domidx;
    domidx.reserve(dims.soc.size() + 2);
    if (dims.nonneg > 0) {
        MSKint64t d;
        check(MSK_appendrplusdomain(task, dims.nonneg, &d), "MSK_appendrplusdomain", task);
        domidx.push_back(d);
    }
    std::vector<std::pair<int, MSKint64t>> soc_domains;
    for (int dim : dims.soc) {
        auto it = std::find_if(soc_domains.begin(), soc_domains.end(), [dim](const auto& e) { return e.first == dim; });
        if (it == soc_domains.end()) {
            MSKint64t d;
            check(MSK_appendquadraticconedomain(task, dim, &d), "MSK_appendquadraticconedomain", task);
            soc_domains.emplace_back(dim, d);
            it = soc_domains.end() - 1;
        }
        domidx.push_back(it->second);
    }
    if (has_quad) {
        MSKint64t d;
        check(MSK_appendrquadraticconedomain(task, m_separable ? 3 : m_r + 2, &d), "MSK_appendrquadraticconedomain", task);
        domidx.insert(domidx.end(), m_separable ? static_cast<std::size_t>(m_r) : 1, d);
    }
    if (!domidx.empty()) check(MSK_appendaccsseq(task, static_cast<MSKint64t>(domidx.size()), domidx.data(), numafe, 0, nullptr), "MSK_appendaccsseq", task);

    m_c.resize(static_cast<std::size_t>(m_n));
    fillValues(P, q, A, b);
    putValues();

    m_xx.resize(nvar);
    m_doty.resize(numafe);
    m_y.resize(m_opts.equality_duals ? m_p : 0);
}

void ConicSolver::update(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b) {
    if (!m_task) throw std::logic_error("mosek::ConicSolver::update() called before setup()");
    if (q.size() != m_n || b.size() != m_p + m_m) throw std::invalid_argument("mosek::ConicSolver::update(): inconsistent dimensions");
    if (!m_native) computeQuadFactor(P, false);
    fillValues(P, q, A, b);
    putValues();
}

// D from P's diagonal, E from the largest |(A D)_ij| per zero-cone row and per cone block.
void ConicSolver::computeScaling(const SparseMatrix<Eigen::ColMajor>& P, const SparseMatrix<Eigen::ColMajor>& A) {
    m_D.setOnes(m_n);
    m_E.setOnes(m_p + m_m);
    if (!m_opts.scaling) return;

    const int* pouter = P.outerIndexPtr();
    const int* pinner = P.innerIndexPtr();
    const double* pval = P.valuePtr();
    for (int j = 0; j < m_n; ++j)
        for (int k = pouter[j]; k < pouter[j + 1]; ++k)
            if (pinner[k] == j && pval[k] > 0) m_D[j] = 1.0 / std::sqrt(pval[k]);

    VectorXr& rowmax = m_E;  // reuse storage: first row maxima, then their inverses
    rowmax.setZero();
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    const double* aval = A.valuePtr();
    for (int j = 0; j < m_n; ++j)
        for (int k = outer[j]; k < outer[j + 1]; ++k) rowmax[inner[k]] = std::max(rowmax[inner[k]], std::abs(aval[k]) * m_D[j]);

    // One common factor per cone block (each nonneg row is its own block).
    int i = m_p + m_dims.nonneg;
    for (int d : m_dims.soc) {
        const double mx = rowmax.segment(i, d).maxCoeff();
        rowmax.segment(i, d).setConstant(mx);
        i += d;
    }
    for (int r = 0; r < rowmax.size(); ++r) rowmax[r] = rowmax[r] > 0 ? 1.0 / rowmax[r] : 1.0;
}

// Copies the (scaled) numerical values into the cached value arrays (in the cached entry order).
void ConicSolver::fillValues(const SparseMatrix<Eigen::ColMajor>& P, const VectorXr& q, const SparseMatrix<Eigen::ColMajor>& A, const VectorXr& b) {
    computeScaling(P, A);
    const double* D = m_D.data();
    const double* E = m_E.data();

    const double* pval = P.valuePtr();
    for (std::size_t e = 0; e < m_q_src.size(); ++e) m_q_val[e] = pval[m_q_src[e]] * D[m_q_subi[e]] * D[m_q_subj[e]];
    for (int j = 0; j < m_n; ++j) m_c[static_cast<std::size_t>(j)] = q[j] * D[j];

    const double* aval = A.valuePtr();
    for (std::size_t e = 0; e < m_lin_src.size(); ++e) m_a_val[e] = aval[m_lin_src[e]] * E[m_a_subi[e]] * D[m_a_subj[e]];
    // F = -A on the cone rows (sign flip of Gx + s = h  ->  -Gx + h in K).
    for (std::size_t e = 0; e < m_cone_src.size(); ++e) m_f_val[e] = -aval[m_cone_src[e]] * E[m_p + m_f_afe[e]] * D[m_f_var[e]];
    std::size_t e = m_cone_src.size();
    if (m_r > 0 && m_separable) {
        for (std::size_t k = 0; k < m_l_val.size(); ++k) {
            m_f_val[e++] = 1.0;  // t_k
            m_f_val[e++] = m_l_val[k] * D[m_l_col[k]];
        }
    } else if (m_r > 0) {
        m_f_val[e++] = 1.0;  // t
        for (std::size_t k = 0; k < m_l_val.size(); ++k) m_f_val[e++] = m_l_val[k] * D[m_l_col[k]];
    }
    for (int r = 0; r < m_p; ++r) m_blc[static_cast<std::size_t>(r)] = b[r] * E[r];
    for (int r = 0; r < m_m; ++r) m_g[static_cast<std::size_t>(r)] = b[m_p + r] * E[m_p + r];  // g = h
}

void ConicSolver::putValues() {
    MSKtask_t task = m_task.get();
    check(MSK_putcslice(task, 0, m_n, m_c.data()), "MSK_putcslice", task);
    if (m_native && !m_q_val.empty())
        check(MSK_putqobj(task, static_cast<MSKint32t>(m_q_val.size()), m_q_subi.data(), m_q_subj.data(), m_q_val.data()), "MSK_putqobj", task);
    if (!m_a_subi.empty()) check(MSK_putaijlist64(task, static_cast<MSKint64t>(m_a_subi.size()), m_a_subi.data(), m_a_subj.data(), m_a_val.data()), "MSK_putaijlist64", task);
    if (m_p > 0) check(MSK_putconboundslice(task, 0, m_p, m_bkc.data(), m_blc.data(), m_blc.data()), "MSK_putconboundslice", task);
    if (!m_f_afe.empty()) check(MSK_putafefentrylist(task, static_cast<MSKint64t>(m_f_afe.size()), m_f_afe.data(), m_f_var.data(), m_f_val.data()), "MSK_putafefentrylist", task);
    if (!m_g.empty()) check(MSK_putafegslice(task, 0, static_cast<MSKint64t>(m_g.size()), m_g.data()), "MSK_putafegslice", task);
}

Status ConicSolver::optimize() {
    if (!m_task) throw std::logic_error("mosek::ConicSolver::optimize() called before setup()");
    MSKtask_t task = m_task.get();

    m_trm = MSK_RES_OK;
    check(MSK_optimizetrm(task, &m_trm), "MSK_optimizetrm", task);

    MSKsolstae solsta = MSK_SOL_STA_UNKNOWN;
    check(MSK_getsolsta(task, MSK_SOL_ITR, &solsta), "MSK_getsolsta", task);

    MSKint32t iters = 0;
    check(MSK_getintinf(task, MSK_IINF_INTPNT_ITER, &iters), "MSK_getintinf(INTPNT_ITER)", task);
    m_iterations = iters;
    check(MSK_getdouinf(task, MSK_DINF_OPTIMIZER_TIME, &m_optimizer_time), "MSK_getdouinf(OPTIMIZER_TIME)", task);

    if (solsta == MSK_SOL_STA_OPTIMAL)
        m_status = Status::Optimal;
    else if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER)
        m_status = Status::PrimalInfeasible;
    else if (solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
        m_status = Status::DualInfeasible;
    else if (m_trm == MSK_RES_TRM_MAX_ITERATIONS)
        m_status = Status::IterationLimit;
    else if (m_trm == MSK_RES_TRM_MAX_TIME)
        m_status = Status::TimeLimit;
    else if (solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS)
        m_status = Status::NearOptimal;
    else if (m_trm == MSK_RES_TRM_STALL || m_trm == MSK_RES_TRM_NUMERICAL_PROBLEM)
        m_status = Status::NumericalFailure;
    else
        m_status = Status::Unknown;
    return m_status;
}

void ConicSolver::extract() {
    MSKtask_t task = m_task.get();
    // Undo the scaling: x = D xh, z = E zh.
    check(MSK_getxx(task, MSK_SOL_ITR, m_xx.data()), "MSK_getxx", task);
    m_xx.head(m_n).array() *= m_D.array();
    if (m_doty.size() > 0) {
        check(MSK_getaccdotys(task, MSK_SOL_ITR, m_doty.data()), "MSK_getaccdotys", task);
        m_doty.head(m_m).array() *= m_E.tail(m_m).array();
    }
    if (m_y.size() > 0) {
        check(MSK_gety(task, MSK_SOL_ITR, m_y.data()), "MSK_gety", task);
        // MOSEK: c - A^T y - F^T doty = 0; Clarabel/QOCO: P x + q + A^T z = 0
        m_y.array() *= -m_E.head(m_p).array();
    }
    check(MSK_getprimalobj(task, MSK_SOL_ITR, &m_objective), "MSK_getprimalobj", task);
}

}  // namespace cardillo::solver::mosek
