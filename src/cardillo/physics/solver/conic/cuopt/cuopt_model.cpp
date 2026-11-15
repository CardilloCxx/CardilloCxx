#include "cuopt_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace cardillo::solver::cuopt {

namespace {

CuoptInt checkedIndex(std::int64_t value, const char* what) {
    if (value > std::numeric_limits<CuoptInt>::max()) throw std::invalid_argument(std::string("cuOpt model: ") + what + " exceeds the index range of cuOpt");
    return static_cast<CuoptInt>(value);
}

bool columnHasNonzero(const CscMatrix& M, Eigen::Index j) {
    for (CscMatrix::InnerIterator it(M, j); it; ++it)
        if (it.value() != 0) return true;
    return false;
}

// Violation of v in the cone of one block: |v| for the zero cone, max(0, -v) for R+ and
// max(0, ||tail|| - head) for a second-order cone.
real_t coneViolation(const Eigen::Ref<const VectorXr>& v, bool zero_cone) {
    if (zero_cone) return v.lpNorm<Eigen::Infinity>();
    const real_t tail = v.size() > 1 ? v.tail(v.size() - 1).norm() : real_t(0);
    return std::max(real_t(0), tail - v[0]);
}

}  // namespace

int ConeDims::rows() const {
    return zero + nonneg + std::accumulate(soc.begin(), soc.end(), 0);
}

void buildModel(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims, Model& model) {
    const Eigen::Index n = P.cols();
    const Eigen::Index m = A.rows();
    if (P.rows() != n || q.size() != n || A.cols() != n || b.size() != m) throw std::invalid_argument("cuOpt model: inconsistent dimensions of P, q, A and b");
    if (dims.zero < 0 || dims.nonneg < 0 || std::any_of(dims.soc.begin(), dims.soc.end(), [](int d) { return d < 1; })) throw std::invalid_argument("cuOpt model: negative or empty cone dimension");
    if (dims.rows() != m) throw std::invalid_argument("cuOpt model: the cone dimensions do not add up to the rows of A");

    // Full symmetric P from its upper triangle, without explicit zeros.
    std::vector<Eigen::Triplet<real_t>> triplets;
    triplets.reserve(2 * static_cast<std::size_t>(P.nonZeros()));
    for (Eigen::Index j = 0; j < n; ++j)
        for (CscMatrix::InnerIterator it(P, j); it; ++it) {
            if (it.row() > j || it.value() == 0) continue;
            triplets.emplace_back(it.row(), j, it.value());
            if (it.row() != j) triplets.emplace_back(j, it.row(), it.value());
        }
    CscMatrix P_full(n, n);
    P_full.setFromTriplets(triplets.begin(), triplets.end());

    // x_j enters the model iff its column of P is nonzero.
    model.x_var.assign(static_cast<std::size_t>(n), -1);
    std::int64_t nx = 0;
    for (Eigen::Index j = 0; j < n; ++j) {
        if (P_full.outerIndexPtr()[j + 1] > P_full.outerIndexPtr()[j]) {
            model.x_var[static_cast<std::size_t>(j)] = static_cast<CuoptInt>(nx++);
            continue;
        }
        if (columnHasNonzero(A, j)) throw std::invalid_argument("cuOpt model: x_" + std::to_string(j) + " has a zero column of P but enters A; the Wolfe dual does not determine it");
        if (q[j] != 0) throw std::invalid_argument("cuOpt model: x_" + std::to_string(j) + " has a linear cost but no curvature and no constraint; the cone program is unbounded");
    }
    model.z_offset = checkedIndex(nx, "the number of variables");
    model.num_vars = checkedIndex(nx + m, "the number of variables");
    model.num_rows = model.z_offset;

    // Row j of [P A'] = column j of P_full followed by column j of A (the rows of A are the z indices).
    model.row_offsets.clear();
    model.columns.clear();
    model.values.clear();
    model.rhs.clear();
    model.row_offsets.reserve(static_cast<std::size_t>(nx) + 1);
    model.columns.reserve(static_cast<std::size_t>(P_full.nonZeros() + A.nonZeros()));
    model.values.reserve(static_cast<std::size_t>(P_full.nonZeros() + A.nonZeros()));
    model.rhs.reserve(static_cast<std::size_t>(nx));
    model.row_offsets.push_back(0);
    for (Eigen::Index j = 0; j < n; ++j) {
        if (model.x_var[static_cast<std::size_t>(j)] < 0) continue;
        for (CscMatrix::InnerIterator it(P_full, j); it; ++it) {
            model.columns.push_back(model.x_var[static_cast<std::size_t>(it.row())]);
            model.values.push_back(it.value());
        }
        for (CscMatrix::InnerIterator it(A, j); it; ++it) {
            if (it.value() == 0) continue;
            model.columns.push_back(model.z_offset + static_cast<CuoptInt>(it.row()));
            model.values.push_back(it.value());
        }
        model.row_offsets.push_back(checkedIndex(static_cast<std::int64_t>(model.columns.size()), "the number of nonzeros"));
        model.rhs.push_back(-q[j]);
    }

    // Quadratic objective Q = P/2 over both triangles of P.
    model.quad_rows.clear();
    model.quad_cols.clear();
    model.quad_values.clear();
    for (Eigen::Index j = 0; j < n; ++j)
        for (CscMatrix::InnerIterator it(P_full, j); it; ++it) {
            model.quad_rows.push_back(model.x_var[static_cast<std::size_t>(it.row())]);
            model.quad_cols.push_back(model.x_var[static_cast<std::size_t>(j)]);
            model.quad_values.push_back(real_t(0.5) * it.value());
        }

    // Bounds: x and the zero-cone multipliers free, R+ multipliers and cone heads >= 0. Linear
    // objective (0, b).
    constexpr real_t inf = std::numeric_limits<real_t>::infinity();
    model.lower.assign(static_cast<std::size_t>(model.num_vars), -inf);
    model.upper.assign(static_cast<std::size_t>(model.num_vars), inf);
    model.objective.assign(static_cast<std::size_t>(model.num_vars), 0.0);
    for (Eigen::Index i = 0; i < m; ++i) model.objective[static_cast<std::size_t>(model.z_offset + i)] = b[i];
    CuoptInt row = model.z_offset + dims.zero;
    for (int i = 0; i < dims.nonneg; ++i) model.lower[static_cast<std::size_t>(row++)] = 0.0;

    model.cone_offsets.assign(1, 0);
    model.cone_vars.clear();
    model.cone_coeffs.clear();
    for (int d : dims.soc) {
        model.lower[static_cast<std::size_t>(row)] = 0.0;
        if (d >= 2) {
            model.cone_vars.push_back(row);
            model.cone_coeffs.push_back(-1.0);
            for (int k = 1; k < d; ++k) {
                model.cone_vars.push_back(row + k);
                model.cone_coeffs.push_back(1.0);
            }
            model.cone_offsets.push_back(static_cast<CuoptInt>(model.cone_vars.size()));
        }
        row += d;
    }
}

void recoverSolution(const Model& model, const CscMatrix& A, const VectorXr& b, const std::vector<real_t>& primal, VectorXr& x, VectorXr& z, VectorXr& s) {
    if (static_cast<std::int64_t>(primal.size()) != model.num_vars) throw std::invalid_argument("cuOpt model: the primal solution does not match the model");
    const Eigen::Index n = static_cast<Eigen::Index>(model.x_var.size());
    x.resize(n);
    for (Eigen::Index j = 0; j < n; ++j) {
        const CuoptInt v = model.x_var[static_cast<std::size_t>(j)];
        x[j] = v >= 0 ? primal[static_cast<std::size_t>(v)] : real_t(0);
    }
    z = Eigen::Map<const VectorXr>(primal.data() + model.z_offset, A.rows());
    s = b;
    s.noalias() -= A * x;
}

real_t KktResiduals::max() const {
    return std::max({primal, dual, dual_cone, complementarity});
}

KktResiduals kktResiduals(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims, const VectorXr& x, const VectorXr& s, const VectorXr& z) {
    KktResiduals r;
    VectorXr stationarity = q;
    stationarity.noalias() += P.selfadjointView<Eigen::Upper>() * x;
    stationarity.noalias() += A.transpose() * z;
    r.dual = stationarity.lpNorm<Eigen::Infinity>() / (1 + q.lpNorm<Eigen::Infinity>());

    real_t s_violation = 0, z_violation = 0, complementarity = 0;
    auto block = [&](Eigen::Index start, Eigen::Index size, bool zero_cone) {
        const auto sb = s.segment(start, size);
        const auto zb = z.segment(start, size);
        s_violation = std::max(s_violation, coneViolation(sb, zero_cone));
        if (zero_cone) return;  // z is free on the zero cone and s'z vanishes with s
        z_violation = std::max(z_violation, coneViolation(zb, false));
        complementarity = std::max(complementarity, std::abs(sb.dot(zb)) / (1 + sb.norm() * zb.norm()));
    };
    Eigen::Index row = 0;
    if (dims.zero > 0) block(row, dims.zero, true);
    row += dims.zero;
    for (int i = 0; i < dims.nonneg; ++i) block(row++, 1, false);
    for (int d : dims.soc) {
        block(row, d, false);
        row += d;
    }
    r.primal = s_violation / (1 + b.lpNorm<Eigen::Infinity>());
    r.dual_cone = z_violation / (1 + z.lpNorm<Eigen::Infinity>());
    r.complementarity = complementarity;
    return r;
}

}  // namespace cardillo::solver::cuopt
