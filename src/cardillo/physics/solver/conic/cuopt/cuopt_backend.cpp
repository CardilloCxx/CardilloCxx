#include "cuopt_backend.hpp"

#include <cuopt/mathematical_optimization/cuopt_c.h>

#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace cardillo::solver::cuopt {

// The arrays of Model are handed to cuOpt without conversion. The constructor checks the same for the
// loaded library (cuOptGetIntSize, cuOptGetFloatSize).
static_assert(std::is_same_v<cuopt_int_t, CuoptInt>, "the cuOpt header must use CuoptInt indices (cuOpt's default build)");
static_assert(std::is_same_v<cuopt_float_t, real_t>, "the cuOpt header must use real_t as floating-point type");

namespace {

const char* returnCodeName(cuopt_int_t code) {
    switch (code) {
        case CUOPT_SUCCESS:
            return "CUOPT_SUCCESS";
        case CUOPT_INVALID_ARGUMENT:
            return "CUOPT_INVALID_ARGUMENT";
        case CUOPT_MPS_FILE_ERROR:
            return "CUOPT_MPS_FILE_ERROR";
        case CUOPT_MPS_PARSE_ERROR:
            return "CUOPT_MPS_PARSE_ERROR";
        case CUOPT_VALIDATION_ERROR:
            return "CUOPT_VALIDATION_ERROR";
        case CUOPT_OUT_OF_MEMORY:
            return "CUOPT_OUT_OF_MEMORY";
        case CUOPT_RUNTIME_ERROR:
            return "CUOPT_RUNTIME_ERROR";
        default:
            return "unknown cuOpt return code";
    }
}

const char* terminationName(cuopt_int_t status) {
    switch (status) {
        case CUOPT_TERMINATION_STATUS_NO_TERMINATION:
            return "NO_TERMINATION";
        case CUOPT_TERMINATION_STATUS_OPTIMAL:
            return "OPTIMAL";
        case CUOPT_TERMINATION_STATUS_INFEASIBLE:
            return "INFEASIBLE";
        case CUOPT_TERMINATION_STATUS_UNBOUNDED:
            return "UNBOUNDED";
        case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
            return "ITERATION_LIMIT";
        case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
            return "TIME_LIMIT";
        case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
            return "NUMERICAL_ERROR";
        case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
            return "PRIMAL_FEASIBLE";
        case CUOPT_TERMINATION_STATUS_FEASIBLE_FOUND:
            return "FEASIBLE_FOUND";
        case CUOPT_TERMINATION_STATUS_CONCURRENT_LIMIT:
            return "CONCURRENT_LIMIT";
        case CUOPT_TERMINATION_STATUS_WORK_LIMIT:
            return "WORK_LIMIT";
        case CUOPT_TERMINATION_STATUS_UNBOUNDED_OR_INFEASIBLE:
            return "UNBOUNDED_OR_INFEASIBLE";
        default:
            return "unknown termination status";
    }
}

Status mapTermination(cuopt_int_t status) {
    switch (status) {
        case CUOPT_TERMINATION_STATUS_OPTIMAL:
            return Status::Optimal;
        case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
        case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
        case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
        case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
            return Status::Inaccurate;
        default:
            return Status::Failed;
    }
}

// The cuOpt problem is the Wolfe dual of the cone program, so its statuses swap meaning.
std::string terminationMessage(cuopt_int_t status) {
    std::string msg = std::string("cuOpt terminated with ") + terminationName(status);
    if (status == CUOPT_TERMINATION_STATUS_UNBOUNDED) msg += " (the dual is unbounded: the cone program is infeasible)";
    if (status == CUOPT_TERMINATION_STATUS_INFEASIBLE) msg += " (the dual is infeasible: the cone program is unbounded or infeasible)";
    return msg;
}

void check(cuopt_int_t code, const char* function) {
    if (code != CUOPT_SUCCESS) throw std::runtime_error(std::string("cuOpt: ") + function + " returned " + returnCodeName(code));
}

cuopt_int_t presolveValue(Presolve p) {
    switch (p) {
        case Presolve::Off:
            return CUOPT_PRESOLVE_OFF;
        case Presolve::Papilo:
            return CUOPT_PRESOLVE_PAPILO;
        case Presolve::Pslp:
            return CUOPT_PRESOLVE_PSLP;
        default:
            return CUOPT_PRESOLVE_DEFAULT;
    }
}

// The destroy functions take the address of the handle and reset it, so the handles are owned by
// small structs instead of std::unique_ptr.
struct SettingsHandle {
    cuOptSolverSettings h{nullptr};
    SettingsHandle() = default;
    SettingsHandle(const SettingsHandle&) = delete;
    SettingsHandle& operator=(const SettingsHandle&) = delete;
    ~SettingsHandle() {
        if (h) cuOptDestroySolverSettings(&h);
    }
};

struct ProblemHandle {
    cuOptOptimizationProblem h{nullptr};
    ProblemHandle() = default;
    ProblemHandle(const ProblemHandle&) = delete;
    ProblemHandle& operator=(const ProblemHandle&) = delete;
    ~ProblemHandle() { reset(); }
    void reset() {
        if (h) cuOptDestroyProblem(&h);
    }
};

struct SolutionHandle {
    cuOptSolution h{nullptr};
    SolutionHandle() = default;
    SolutionHandle(const SolutionHandle&) = delete;
    SolutionHandle& operator=(const SolutionHandle&) = delete;
    ~SolutionHandle() {
        if (h) cuOptDestroySolution(&h);
    }
};

}  // namespace

const char* toString(Status status) {
    switch (status) {
        case Status::Optimal:
            return "optimal";
        case Status::Inaccurate:
            return "inaccurate";
        default:
            return "failed";
    }
}

struct Solver::Impl {
    SettingsHandle settings;
    ProblemHandle problem;
    Model model;
    std::vector<char> senses;
    std::vector<char> types;
    std::vector<real_t> primal;
    VectorXr x, z, s;
    Info info;
};

std::string Solver::version() {
    cuopt_int_t major = 0, minor = 0, patch = 0;
    check(cuOptGetVersion(&major, &minor, &patch), "cuOptGetVersion");
    return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch);
}

Solver::Solver(const Settings& settings) : m_impl(std::make_unique<Impl>()) {
    // Cones arrived in 26.06; the SOCP path received correctness fixes up to the 26.10.1 hotfix.
    cuopt_int_t major = 0, minor = 0, patch = 0;
    check(cuOptGetVersion(&major, &minor, &patch), "cuOptGetVersion");
    if (std::tie(major, minor, patch) < std::make_tuple(26, 10, 1)) throw std::runtime_error("cuOpt " + version() + " is too old; the cuOpt backend needs 26.10.1 or newer");
    if (cuOptGetIntSize() != sizeof(CuoptInt) || cuOptGetFloatSize() != sizeof(real_t))
        throw std::runtime_error("the loaded cuOpt library uses " + std::to_string(cuOptGetIntSize()) + "-byte indices and " + std::to_string(cuOptGetFloatSize()) +
                                 "-byte floats; the cuOpt backend needs " + std::to_string(sizeof(CuoptInt)) + " and " + std::to_string(sizeof(real_t)));

    cuOptSolverSettings h = nullptr;
    check(cuOptCreateSolverSettings(&h), "cuOptCreateSolverSettings");
    m_impl->settings.h = h;
    auto setInt = [h](const char* name, cuopt_int_t value) { check(cuOptSetIntegerParameter(h, name, value), name); };
    auto setFloat = [h](const char* name, cuopt_float_t value) { check(cuOptSetFloatParameter(h, name, value), name); };
    // Quadratic constraints always run the barrier method; set it explicitly.
    setInt(CUOPT_METHOD, CUOPT_METHOD_BARRIER);
    setInt(CUOPT_LOG_TO_CONSOLE, settings.verbose ? 1 : 0);
    setFloat(CUOPT_ABSOLUTE_PRIMAL_TOLERANCE, settings.tol_abs);
    setFloat(CUOPT_ABSOLUTE_DUAL_TOLERANCE, settings.tol_abs);
    setFloat(CUOPT_ABSOLUTE_GAP_TOLERANCE, settings.tol_abs);
    setFloat(CUOPT_RELATIVE_PRIMAL_TOLERANCE, settings.tol_rel);
    setFloat(CUOPT_RELATIVE_DUAL_TOLERANCE, settings.tol_rel);
    setFloat(CUOPT_RELATIVE_GAP_TOLERANCE, settings.tol_rel);
    if (settings.time_limit > 0) setFloat(CUOPT_TIME_LIMIT, settings.time_limit);
    if (settings.iteration_limit >= 0) setInt(CUOPT_ITERATION_LIMIT, settings.iteration_limit);
    setInt(CUOPT_PRESOLVE, presolveValue(settings.presolve));
    setInt(CUOPT_CUDSS_DETERMINISTIC, settings.deterministic ? 1 : 0);
}

Solver::~Solver() = default;

void Solver::setup(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims) {
    Impl& im = *m_impl;
    im.problem.reset();
    buildModel(P, q, A, b, dims, im.model);
    const Model& md = im.model;
    im.senses.assign(static_cast<std::size_t>(md.num_rows), CUOPT_EQUAL);
    im.types.assign(static_cast<std::size_t>(md.num_vars), CUOPT_CONTINUOUS);

    check(cuOptCreateProblem(md.num_rows, md.num_vars, CUOPT_MINIMIZE, 0.0, md.objective.data(), md.row_offsets.data(), md.columns.data(), md.values.data(), im.senses.data(), md.rhs.data(),
                             md.lower.data(), md.upper.data(), im.types.data(), &im.problem.h),
          "cuOptCreateProblem");
    if (!md.quad_values.empty())
        check(cuOptSetQuadraticObjective(im.problem.h, static_cast<cuopt_int_t>(md.quad_values.size()), md.quad_rows.data(), md.quad_cols.data(), md.quad_values.data()), "cuOptSetQuadraticObjective");
    for (int k = 0; k < md.numCones(); ++k) {
        const CuoptInt begin = md.cone_offsets[static_cast<std::size_t>(k)];
        const CuoptInt len = md.cone_offsets[static_cast<std::size_t>(k) + 1] - begin;
        const CuoptInt* vars = md.cone_vars.data() + begin;
        check(cuOptAddQuadraticConstraint(im.problem.h, len, vars, vars, md.cone_coeffs.data() + begin, 0, nullptr, nullptr, CUOPT_LESS_THAN, 0.0), "cuOptAddQuadraticConstraint");
    }
}

Status Solver::solve() {
    Impl& im = *m_impl;
    if (!im.problem.h) throw std::logic_error("cuOpt: solve() without setup()");
    im.info = Info{};

    SolutionHandle sol;
    const cuopt_int_t rc = cuOptSolve(im.problem.h, im.settings.h, &sol.h);
    if (!sol.h) {
        im.info.message = std::string("cuOptSolve returned ") + returnCodeName(rc) + " without a solution";
        if (rc == CUOPT_VALIDATION_ERROR) im.info.message += " (cuOpt reports this also when no usable NVIDIA GPU and driver are found)";
        return im.info.status;
    }
    if (rc != CUOPT_SUCCESS) {
        char text[1024] = {0};
        cuOptGetErrorString(sol.h, text, static_cast<cuopt_int_t>(sizeof text));
        im.info.message = std::string("cuOptSolve returned ") + returnCodeName(rc) + ": " + text;
        return im.info.status;
    }

    cuopt_int_t termination = -1;
    check(cuOptGetTerminationStatus(sol.h, &termination), "cuOptGetTerminationStatus");
    im.info.termination = termination;
    im.info.status = mapTermination(termination);
    cuopt_float_t value = 0;
    if (cuOptGetSolveTime(sol.h, &value) == CUOPT_SUCCESS) im.info.solve_time = value;
    if (cuOptGetObjectiveValue(sol.h, &value) == CUOPT_SUCCESS) im.info.objective = value;
    cuopt_int_t iterations = -1;
    if (cuOptGetSolutionIntAttribute(sol.h, CUOPT_SOLUTION_ATTR_LP_NUM_ITERATIONS, &iterations) == CUOPT_SUCCESS) im.info.iterations = iterations;

    if (im.info.status == Status::Failed) {
        im.info.message = terminationMessage(termination);
        return im.info.status;
    }
    if (im.info.status == Status::Inaccurate) im.info.message = terminationMessage(termination);
    im.primal.resize(static_cast<std::size_t>(im.model.num_vars));
    check(cuOptGetPrimalSolution(sol.h, im.primal.data()), "cuOptGetPrimalSolution");
    return im.info.status;
}

void Solver::recover(const CscMatrix& A, const VectorXr& b) {
    Impl& im = *m_impl;
    if (im.info.status == Status::Failed) throw std::logic_error("cuOpt: recover() after a failed solve");
    recoverSolution(im.model, A, b, im.primal, im.x, im.z, im.s);
}

const VectorXr& Solver::x() const {
    return m_impl->x;
}
const VectorXr& Solver::z() const {
    return m_impl->z;
}
const VectorXr& Solver::s() const {
    return m_impl->s;
}
const Info& Solver::info() const {
    return m_impl->info;
}

bool Solver::writeProblem(const std::string& path) const {
    return m_impl->problem.h && cuOptWriteProblem(m_impl->problem.h, path.c_str(), CUOPT_FILE_FORMAT_MPS) == CUOPT_SUCCESS;
}

bool deviceAvailable(std::string* reason) {
    try {
        // Projection of p = (0, 1, 0) onto L^3: min 1/2 ||x||^2 - p'x s.t. -x + s = 0, s in L^3.
        CscMatrix I(3, 3);
        I.setIdentity();
        const CscMatrix A = -I;
        const VectorXr q = -Vector3r(0, 1, 0);
        const VectorXr b = VectorXr::Zero(3);
        ConeDims dims;
        dims.soc = {3};
        Solver solver{Settings{}};
        solver.setup(I, q, A, b, dims);
        if (solver.solve() != Status::Failed) return true;
        if (reason) *reason = solver.info().message;
    } catch (const std::exception& e) {
        if (reason) *reason = e.what();
    }
    return false;
}

}  // namespace cardillo::solver::cuopt
