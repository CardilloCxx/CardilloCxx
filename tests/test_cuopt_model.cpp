// Translation of cone programs into cuOpt's form (cuopt_model.hpp: Wolfe dual, cones as quadratic
// constraints), checked without cuOpt and without a GPU:
//  1. structure of the model: valid CSR rows, bounds per cone block, one quadratic constraint per
//     cone of dimension >= 2 with coefficients -1 (head) and +1 (tail), left-out variables;
//  2. reference injection: a known KKT point (x, z) of the cone program satisfies the rows, bounds
//     and cones of the model, and its objective is minus the optimal value of the cone program;
//  3. the quadratic objective without factor 1/2: x'Qx = 1/2 x'Px;
//  4. recovery of x, z and s = b - Ax from a primal vector of the model;
//  5. P given as upper triangle or full matrix gives the same model; invalid input throws.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cuopt_test_problems.hpp"

using namespace cuopt_test;
using cardillo::solver::cuopt::buildModel;
using cardillo::solver::cuopt::CuoptInt;
using cardillo::solver::cuopt::kktResiduals;
using cardillo::solver::cuopt::Model;
using cardillo::solver::cuopt::recoverSolution;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) {
        ++g_failures;
        std::cerr << "FAILED: " << what << std::endl;
    }
}

constexpr real_t kTol = 1e-12;

std::vector<real_t> modelPoint(const Model& md, const VectorXr& x, const VectorXr& z) {
    std::vector<real_t> y(static_cast<std::size_t>(md.num_vars), 0.0);
    for (std::size_t j = 0; j < md.x_var.size(); ++j)
        if (md.x_var[j] >= 0) y[static_cast<std::size_t>(md.x_var[j])] = x[static_cast<Eigen::Index>(j)];
    for (Eigen::Index i = 0; i < z.size(); ++i) y[static_cast<std::size_t>(md.z_offset + i)] = z[i];
    return y;
}

real_t maxRowResidual(const Model& md, const std::vector<real_t>& y) {
    real_t r = 0;
    for (int row = 0; row < md.num_rows; ++row) {
        real_t sum = -md.rhs[static_cast<std::size_t>(row)];
        for (auto k = md.row_offsets[static_cast<std::size_t>(row)]; k < md.row_offsets[static_cast<std::size_t>(row) + 1]; ++k)
            sum += md.values[static_cast<std::size_t>(k)] * y[static_cast<std::size_t>(md.columns[static_cast<std::size_t>(k)])];
        r = std::max(r, std::abs(sum));
    }
    return r;
}

// Largest value of the cone constraints sum coeff y^2 (must be <= 0).
real_t maxConeValue(const Model& md, const std::vector<real_t>& y) {
    real_t r = -INFINITY;
    for (int k = 0; k < md.numCones(); ++k) {
        real_t sum = 0;
        for (auto i = md.cone_offsets[static_cast<std::size_t>(k)]; i < md.cone_offsets[static_cast<std::size_t>(k) + 1]; ++i) {
            const real_t v = y[static_cast<std::size_t>(md.cone_vars[static_cast<std::size_t>(i)])];
            sum += md.cone_coeffs[static_cast<std::size_t>(i)] * v * v;
        }
        r = std::max(r, sum);
    }
    return r;
}

real_t quadraticObjective(const Model& md, const std::vector<real_t>& y) {
    real_t sum = 0;
    for (std::size_t k = 0; k < md.quad_values.size(); ++k) sum += md.quad_values[k] * y[static_cast<std::size_t>(md.quad_rows[k])] * y[static_cast<std::size_t>(md.quad_cols[k])];
    return sum;
}

real_t linearObjective(const Model& md, const std::vector<real_t>& y) {
    real_t sum = 0;
    for (std::size_t k = 0; k < y.size(); ++k) sum += md.objective[k] * y[k];
    return sum;
}

bool boundsHold(const Model& md, const std::vector<real_t>& y) {
    for (std::size_t k = 0; k < y.size(); ++k)
        if (y[k] < md.lower[k] || y[k] > md.upper[k]) return false;
    return true;
}

void checkStructure(const Model& md, const KnownProblem& pr, const std::string& tag) {
    const int m = pr.dims.rows();
    int coupled = 0;
    for (std::size_t j = 0; j < md.x_var.size(); ++j) {
        const bool zero_column = pr.P_full.col(static_cast<Eigen::Index>(j)).cwiseAbs().sum() == 0;
        check((md.x_var[j] < 0) == zero_column, tag + ": x_" + std::to_string(j) + " left out iff its column of P is zero");
        if (md.x_var[j] >= 0) check(md.x_var[j] == coupled++, tag + ": x variables numbered consecutively");
    }
    check(md.z_offset == coupled && md.num_vars == coupled + m && md.num_rows == coupled, tag + ": dimensions");
    check(md.row_offsets.size() == static_cast<std::size_t>(md.num_rows) + 1 && md.row_offsets.front() == 0, tag + ": row offsets");
    check(md.row_offsets.back() == static_cast<CuoptInt>(md.columns.size()) && md.columns.size() == md.values.size() && md.rhs.size() == static_cast<std::size_t>(md.num_rows),
          tag + ": CSR array sizes");
    bool sorted = true, in_range = true, nonzero = true;
    for (int row = 0; row < md.num_rows; ++row)
        for (auto k = md.row_offsets[static_cast<std::size_t>(row)]; k < md.row_offsets[static_cast<std::size_t>(row) + 1]; ++k) {
            const auto c = md.columns[static_cast<std::size_t>(k)];
            in_range &= c >= 0 && c < md.num_vars;
            nonzero &= md.values[static_cast<std::size_t>(k)] != 0;
            if (k > md.row_offsets[static_cast<std::size_t>(row)]) sorted &= c > md.columns[static_cast<std::size_t>(k) - 1];
        }
    check(sorted && in_range, tag + ": CSR column indices sorted, unique and in range");
    check(nonzero, tag + ": no explicit zeros in the rows");

    const std::size_t nv = static_cast<std::size_t>(md.num_vars);
    check(md.lower.size() == nv && md.upper.size() == nv && md.objective.size() == nv, tag + ": bound and objective sizes");
    bool bounds = std::all_of(md.upper.begin(), md.upper.end(), [](real_t u) { return std::isinf(u) && u > 0; });
    for (int k = 0; k < md.z_offset; ++k) bounds &= std::isinf(md.lower[static_cast<std::size_t>(k)]);
    int row = md.z_offset;
    for (int i = 0; i < pr.dims.zero; ++i) bounds &= std::isinf(md.lower[static_cast<std::size_t>(row++)]);
    for (int i = 0; i < pr.dims.nonneg; ++i) bounds &= md.lower[static_cast<std::size_t>(row++)] == 0;
    int cones = 0;
    bool cone_rows = true;
    for (int d : pr.dims.soc) {
        bounds &= md.lower[static_cast<std::size_t>(row)] == 0;
        for (int i = 1; i < d; ++i) bounds &= std::isinf(md.lower[static_cast<std::size_t>(row + i)]);
        if (d >= 2) {
            const auto begin = md.cone_offsets[static_cast<std::size_t>(cones)];
            cone_rows &= md.cone_offsets[static_cast<std::size_t>(cones) + 1] - begin == d;
            for (int i = 0; i < d && cone_rows; ++i) {
                cone_rows &= md.cone_vars[static_cast<std::size_t>(begin + i)] == row + i;
                cone_rows &= md.cone_coeffs[static_cast<std::size_t>(begin + i)] == (i == 0 ? -1.0 : 1.0);
            }
            ++cones;
        }
        row += d;
    }
    check(bounds, tag + ": free x and zero-cone multipliers, R+ multipliers and cone heads >= 0, tails free");
    check(md.numCones() == cones && cone_rows, tag + ": one quadratic constraint per cone of dimension >= 2, head -1, tail +1");

    bool quad_in_x = true;
    for (std::size_t k = 0; k < md.quad_values.size(); ++k) quad_in_x &= md.quad_rows[k] < md.z_offset && md.quad_cols[k] < md.z_offset;
    int nnz = 0;
    for (Eigen::Index j = 0; j < pr.P_full.cols(); ++j)
        for (CscMatrix::InnerIterator it(pr.P_full, j); it; ++it) nnz += it.value() != 0;
    check(quad_in_x && static_cast<int>(md.quad_values.size()) == nnz, tag + ": quadratic objective = nonzeros of P on the x block, each once");
}

// The known KKT point injected into the model.
void checkInjection(const Model& md, const KnownProblem& pr, const std::string& tag) {
    const std::vector<real_t> y = modelPoint(md, pr.x, pr.z);
    const real_t scale = 1 + pr.q.lpNorm<Eigen::Infinity>() + pr.z.lpNorm<Eigen::Infinity>();
    check(maxRowResidual(md, y) < kTol * scale, tag + ": KKT point satisfies the rows P x + A'z = -q");
    check(boundsHold(md, y), tag + ": KKT point satisfies the bounds");
    if (md.numCones() > 0) check(maxConeValue(md, y) < kTol * scale * scale, tag + ": KKT point satisfies the cone constraints");
    const real_t quad = quadraticObjective(md, y);
    const real_t half = 0.5 * pr.x.dot(pr.P_full * pr.x);
    check(std::abs(quad - half) < kTol * (1 + std::abs(half)), tag + ": x'Qx = 1/2 x'Px (Q = P/2)");
    const real_t dual = quad + linearObjective(md, y);
    check(std::abs(dual + pr.value) < kTol * (1 + std::abs(pr.value)), tag + ": dual objective " + std::to_string(dual) + " = -optimal value " + std::to_string(-pr.value));

    VectorXr x, z, s;
    recoverSolution(md, pr.A, pr.b, y, x, z, s);
    check((x - pr.x).lpNorm<Eigen::Infinity>() < kTol && (z - pr.z).lpNorm<Eigen::Infinity>() < kTol && (s - pr.s).lpNorm<Eigen::Infinity>() < kTol * scale, tag + ": recovery of the KKT point");
    const auto r = kktResiduals(pr.P, pr.q, pr.A, pr.b, pr.dims, x, s, z);
    check(r.max() < kTol * scale, tag + ": KKT residuals of the KKT point vanish (" + std::to_string(r.max()) + ")");
    if (z.size() > 0) {
        VectorXr z_bad = z;
        z_bad[z.size() - 1] += 0.1;
        check(kktResiduals(pr.P, pr.q, pr.A, pr.b, pr.dims, x, s, z_bad).max() > 1e-3, tag + ": KKT residuals detect a perturbed z");
    }
}

void checkProblem(const KnownProblem& pr, const std::string& tag) {
    Model md;
    buildModel(pr.P, pr.q, pr.A, pr.b, pr.dims, md);
    checkStructure(md, pr, tag);
    checkInjection(md, pr, tag);

    // Recovery places every model variable: y_k = k + 1/2.
    std::vector<real_t> y(static_cast<std::size_t>(md.num_vars));
    for (std::size_t k = 0; k < y.size(); ++k) y[k] = static_cast<real_t>(k) + 0.5;
    VectorXr x, z, s;
    recoverSolution(md, pr.A, pr.b, y, x, z, s);
    bool placed = true;
    for (std::size_t j = 0; j < md.x_var.size(); ++j) placed &= x[static_cast<Eigen::Index>(j)] == (md.x_var[j] >= 0 ? md.x_var[j] + 0.5 : 0.0);
    for (Eigen::Index i = 0; i < z.size(); ++i) placed &= z[i] == md.z_offset + i + 0.5;
    check(placed && (s - (pr.b - pr.A * x)).lpNorm<Eigen::Infinity>() < kTol * (1 + s.lpNorm<Eigen::Infinity>()), tag + ": recovery places x, z and s = b - Ax");

    // P as full matrix: the same model.
    Model full;
    buildModel(pr.P_full, pr.q, pr.A, pr.b, pr.dims, full);
    check(full.row_offsets == md.row_offsets && full.columns == md.columns && full.values == md.values && full.quad_values == md.quad_values && full.quad_rows == md.quad_rows &&
              full.quad_cols == md.quad_cols && full.x_var == md.x_var,
          tag + ": P as upper triangle or full matrix gives the same model");
}

template <class F>
void checkThrows(F f, const std::string& what) {
    try {
        f();
    } catch (const std::invalid_argument&) {
        return;
    }
    check(false, what);
}

void checkInvalidInput() {
    ConeDims dims;
    dims.nonneg = 1;
    dims.soc = {3};
    const KnownProblem pr = randomProblem(7, 4, 1, dims, true);
    Model md;
    ConeDims wrong = dims;
    wrong.soc = {2};
    checkThrows([&] { buildModel(pr.P, pr.q, pr.A, pr.b, wrong, md); }, "cone dimensions that do not match A throw");
    VectorXr q = pr.q;
    q[4] = 1;  // the decoupled variable gets a cost
    checkThrows([&] { buildModel(pr.P, q, pr.A, pr.b, dims, md); }, "a cost on a variable without curvature and constraint throws");
    CscMatrix A = pr.A;
    A.coeffRef(0, 4) = 1;  // the decoupled variable enters A
    checkThrows([&] { buildModel(pr.P, pr.q, A, pr.b, dims, md); }, "a variable without curvature in A throws");
}

}  // namespace

int main() {
    try {
        ConeDims mixed;
        mixed.zero = 3;
        mixed.nonneg = 4;
        mixed.soc = {3, 3, 3, 3, 2, 1, 3, 2, 1};
        checkProblem(randomProblem(1, 12, 2, mixed, true), "mixed, diagonal P");
        checkProblem(randomProblem(2, 12, 2, mixed, false), "mixed, coupled P");

        ConeDims contacts;
        contacts.soc.assign(8, 3);
        checkProblem(randomProblem(3, 18, 0, contacts, true), "frictional contacts only");

        ConeDims bilateral;
        bilateral.zero = 5;
        checkProblem(randomProblem(4, 9, 3, bilateral, true), "bilateral rows only");

        // Closed-form projections onto L^3, L^2 and L^1: inside, polar cone (apex) and boundary.
        for (const VectorXr& p : {VectorXr(Vector3r(2, 1, 0)), VectorXr(Vector3r(-2, 1, 0)), VectorXr(Vector3r(0, 1, 0)), VectorXr(Eigen::Vector2d(0, 1)), VectorXr(Eigen::Vector2d(-1, 3)),
                                  VectorXr(VectorXr::Constant(1, -1)), VectorXr(VectorXr::Constant(1, 2))}) {
            const KnownProblem pr = coneProjection(p);
            std::string tag = "projection of (";
            for (Eigen::Index i = 0; i < p.size(); ++i) tag += (i ? ", " : "") + std::to_string(p[i]);
            checkProblem(pr, tag + ")");
        }
        const KnownProblem boundary = coneProjection(Vector3r(0, 1, 0));
        check((boundary.x - Vector3r(0.5, 0.5, 0)).norm() < kTol && (boundary.z - Vector3r(0.5, -0.5, 0)).norm() < kTol, "closed-form projection of (0, 1, 0)");

        checkInvalidInput();
    } catch (const std::exception& ex) {
        check(false, std::string("exception: ") + ex.what());
    }
    if (g_failures == 0) std::cout << "test_cuopt_model: all checks passed" << std::endl;
    return g_failures == 0 ? 0 : 1;
}
