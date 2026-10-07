#include "scs_backend.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <type_traits>

#include <scs.h>

namespace cardillo::solver::scs {

// The problem data is handed to SCS without copies or conversions, so SCS must be built with its
// default index/scalar types (DLONG and SFLOAT off, enforced in src/CMakeLists.txt).
static_assert(std::is_same_v<scs_int, SpMat::StorageIndex>, "scs_int must match the Eigen sparse index type (build SCS with DLONG=OFF)");
static_assert(std::is_same_v<scs_float, real_t>, "scs_float must match real_t (build SCS with SFLOAT=OFF)");

int ConeDims::rows() const {
    int r = zero + nonneg;
    for (int k : soc) r += k;
    return r;
}

const char* toString(Status s) {
    switch (s) {
        case Status::Solved:
            return "Solved";
        case Status::SolvedInaccurate:
            return "Solved Inaccurate";
        case Status::Infeasible:
            return "Infeasible";
        case Status::Unbounded:
            return "Unbounded";
        case Status::Failed:
            return "Failed";
        case Status::Interrupted:
            return "Interrupted";
        default:
            return "Unknown";
    }
}

namespace {

Status fromScs(scs_int v) {
    switch (v) {
        case SCS_SOLVED:
            return Status::Solved;
        case SCS_SOLVED_INACCURATE:
            return Status::SolvedInaccurate;
        case SCS_INFEASIBLE:
        case SCS_INFEASIBLE_INACCURATE:
            return Status::Infeasible;
        case SCS_UNBOUNDED:
        case SCS_UNBOUNDED_INACCURATE:
            return Status::Unbounded;
        case SCS_SIGINT:
            return Status::Interrupted;
        case SCS_FAILED:
        case SCS_INDETERMINATE:
            return Status::Failed;
        default:
            return Status::Unknown;
    }
}

ScsMatrix view(const SpMat& M) {
    ScsMatrix out;
    out.x = const_cast<scs_float*>(M.valuePtr());
    out.i = const_cast<scs_int*>(M.innerIndexPtr());
    out.p = const_cast<scs_int*>(M.outerIndexPtr());
    out.m = static_cast<scs_int>(M.rows());
    out.n = static_cast<scs_int>(M.cols());
    return out;
}

}  // namespace

Solver::Solver() = default;

Solver::~Solver() { release(); }

void Solver::release() {
    if (m_work) {
        scs_finish(m_work);
        m_work = nullptr;
    }
}

bool Solver::setup(const SpMat& P, const VectorXr& q, const SpMat& A, const VectorXr& b, const ConeDims& dims, const Settings& settings) {
    release();

    const int n = static_cast<int>(P.cols());
    const int m = static_cast<int>(A.rows());
    if (A.cols() != n || q.size() != n || b.size() != m || dims.rows() != m) return false;

    SpMat Pu = P.triangularView<Eigen::Upper>();
    Pu.makeCompressed();
    SpMat Ac;
    const SpMat* Ap = &A;
    if (!A.isCompressed()) {
        Ac = A;
        Ac.makeCompressed();
        Ap = &Ac;
    }

    ScsMatrix Pm = view(Pu);
    ScsMatrix Am = view(*Ap);
    VectorXr qc = q, bc = b;  // SCS takes non-const pointers (it copies them anyway)

    ScsData data;
    data.m = m;
    data.n = n;
    data.A = &Am;
    data.P = Pu.nonZeros() > 0 ? &Pm : nullptr;
    data.b = bc.data();
    data.c = qc.data();

    std::vector<scs_int> soc(dims.soc.begin(), dims.soc.end());
    ScsCone cone{};  // zero-initializes every cone type not used here (box, PSD, exp, power)
    cone.z = dims.zero;
    cone.l = dims.nonneg;
    cone.q = soc.empty() ? nullptr : soc.data();
    cone.qsize = static_cast<scs_int>(soc.size());

    ScsSettings st;
    scs_set_default_settings(&st);
    if (settings.eps_abs >= 0) st.eps_abs = settings.eps_abs;
    if (settings.eps_rel >= 0) st.eps_rel = settings.eps_rel;
    if (settings.eps_infeas >= 0) st.eps_infeas = settings.eps_infeas;
    if (settings.max_iters > 0) st.max_iters = settings.max_iters;
    if (settings.time_limit_secs >= 0) st.time_limit_secs = settings.time_limit_secs;
    if (settings.alpha > 0) st.alpha = settings.alpha;
    if (settings.scale > 0) st.scale = settings.scale;
    if (settings.rho_x > 0) st.rho_x = settings.rho_x;
    if (settings.acceleration_lookback >= 0) st.acceleration_lookback = settings.acceleration_lookback;
    if (settings.acceleration_interval > 0) st.acceleration_interval = settings.acceleration_interval;
    st.adaptive_scale = settings.adaptive_scale ? 1 : 0;
    st.normalize = settings.normalize ? 1 : 0;
    st.verbose = settings.verbose ? 1 : 0;

    m_work = scs_init(&data, &cone, &st);
    if (!m_work) return false;

    m_x.setZero(n);
    m_y.setZero(m);
    m_s.setZero(m);
    m_info = Info{};
    return true;
}

Status Solver::solve(bool warm_start) {
    if (!m_work) return Status::Failed;

    ScsSolution sol;
    sol.x = m_x.data();
    sol.y = m_y.data();
    sol.s = m_s.data();
    ScsInfo info{};
    const scs_int flag = scs_solve(m_work, &sol, &info, warm_start ? 1 : 0);

    m_info.status = fromScs(flag);
    m_info.iterations = static_cast<int>(info.iter);
    m_info.setup_ms = info.setup_time;
    m_info.solve_ms = info.solve_time;
    m_info.lin_sys_ms = info.lin_sys_time;
    m_info.cone_ms = info.cone_time;
    m_info.accel_ms = info.accel_time;
    m_info.accel_accepted = static_cast<int>(info.accepted_accel_steps);
    m_info.accel_rejected = static_cast<int>(info.rejected_accel_steps);
    m_info.scale_updates = static_cast<int>(info.scale_updates);
    m_info.scale = info.scale;
    m_info.res_pri = info.res_pri;
    m_info.res_dual = info.res_dual;
    m_info.gap = info.gap;
    m_info.pobj = info.pobj;
    return m_info.status;
}

namespace {

constexpr char kDumpMagic[8] = {'C', 'C', 'X', 'C', 'O', 'N', 'E', '1'};

template <class T>
void put(std::ofstream& f, const T* p, std::size_t n) {
    f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n * sizeof(T)));
}
template <class T>
void put(std::ofstream& f, T v) {
    put(f, &v, 1);
}
template <class T>
bool get(std::ifstream& f, T* p, std::size_t n) {
    return static_cast<bool>(f.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(n * sizeof(T))));
}
template <class T>
bool get(std::ifstream& f, T& v) {
    return get(f, &v, 1);
}

void putMatrix(std::ofstream& f, const SpMat& M) {
    SpMat C = M;
    C.makeCompressed();
    put<std::int32_t>(f, static_cast<std::int32_t>(C.nonZeros()));
    put(f, C.outerIndexPtr(), static_cast<std::size_t>(C.cols()) + 1);
    put(f, C.innerIndexPtr(), static_cast<std::size_t>(C.nonZeros()));
    put(f, C.valuePtr(), static_cast<std::size_t>(C.nonZeros()));
}

bool getMatrix(std::ifstream& f, SpMat& M, int rows, int cols) {
    std::int32_t nnz = 0;
    if (!get(f, nnz) || nnz < 0) return false;
    std::vector<int> outer(static_cast<std::size_t>(cols) + 1), inner(static_cast<std::size_t>(nnz));
    std::vector<real_t> val(static_cast<std::size_t>(nnz));
    if (!get(f, outer.data(), outer.size()) || !get(f, inner.data(), inner.size()) || !get(f, val.data(), val.size())) return false;
    M = Eigen::Map<const SpMat>(rows, cols, nnz, outer.data(), inner.data(), val.data());
    return true;
}

}  // namespace

bool writeDump(const std::string& path, const ProblemDump& d) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const int n = static_cast<int>(d.P.cols());
    const int m = static_cast<int>(d.A.rows());
    put(f, kDumpMagic, sizeof kDumpMagic);
    for (int v : {n, m, d.dims.zero, d.dims.nonneg, static_cast<int>(d.dims.soc.size())}) put<std::int32_t>(f, v);
    put(f, d.dims.soc.data(), d.dims.soc.size());
    putMatrix(f, d.P);
    putMatrix(f, d.A);
    put(f, d.q.data(), static_cast<std::size_t>(n));
    put(f, d.b.data(), static_cast<std::size_t>(m));
    put<std::int32_t>(f, d.warm ? 1 : 0);
    put<double>(f, d.scale);
    if (d.warm) {
        put(f, d.x0.data(), static_cast<std::size_t>(n));
        put(f, d.y0.data(), static_cast<std::size_t>(m));
        put(f, d.s0.data(), static_cast<std::size_t>(m));
    }
    return static_cast<bool>(f);
}

bool readDump(const std::string& path, ProblemDump& d) {
    std::ifstream f(path, std::ios::binary);
    char magic[8];
    if (!f || !get(f, magic, 8) || std::memcmp(magic, kDumpMagic, 8) != 0) return false;
    std::int32_t hdr[5];
    if (!get(f, hdr, 5)) return false;
    const int n = hdr[0], m = hdr[1];
    d.dims.zero = hdr[2];
    d.dims.nonneg = hdr[3];
    d.dims.soc.resize(static_cast<std::size_t>(hdr[4]));
    if (!get(f, d.dims.soc.data(), d.dims.soc.size())) return false;
    if (!getMatrix(f, d.P, n, n) || !getMatrix(f, d.A, m, n)) return false;
    d.q.resize(n);
    d.b.resize(m);
    if (!get(f, d.q.data(), static_cast<std::size_t>(n)) || !get(f, d.b.data(), static_cast<std::size_t>(m))) return false;
    std::int32_t warm = 0;
    if (!get(f, warm) || !get(f, d.scale)) return false;
    d.warm = warm != 0;
    if (d.warm) {
        d.x0.resize(n);
        d.y0.resize(m);
        d.s0.resize(m);
        if (!get(f, d.x0.data(), static_cast<std::size_t>(n)) || !get(f, d.y0.data(), static_cast<std::size_t>(m)) || !get(f, d.s0.data(), static_cast<std::size_t>(m))) return false;
    }
    return d.dims.rows() == m;
}

void projectOntoCones(VectorXr& v, const ConeDims& dims) {
    int i = dims.zero;
    for (int k = 0; k < dims.nonneg; ++k, ++i) v[i] = std::max(v[i], (real_t)0);
    for (int k : dims.soc) {
        const real_t t = v[i];
        auto x = v.segment(i + 1, k - 1);
        const real_t nx = x.norm();
        if (nx <= t) {
            // inside the cone
        } else if (nx <= -t) {
            v.segment(i, k).setZero();  // inside the polar cone
        } else {
            const real_t a = (real_t)0.5 * (t + nx);
            v[i] = a;
            x *= a / nx;
        }
        i += k;
    }
}

}  // namespace cardillo::solver::scs
