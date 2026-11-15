#pragma once

#include <cstdint>
#include <vector>

#include "../../../../misc/types.hpp"

namespace cardillo::solver::cuopt {

/// Index type of the arrays handed to cuOpt: cuopt_int_t of cuOpt's default build. It is independent of
/// the index type of the Eigen matrices (the model renumbers all indices); cuopt_backend.cpp checks it
/// against the cuOpt header and library.
using CuoptInt = std::int32_t;

/// Cones of the rows of a cone program, in this order: zero cone, nonnegative orthant, then
/// second-order cones L^d = {(t, u) : ||u|| <= t} (the order of ConeProgramAssembler and Clarabel).
struct ConeDims {
    int zero{0};
    int nonneg{0};
    std::vector<int> soc;

    int rows() const;
};

/**
 * @brief The cone program
 *     min 1/2 x'Px + q'x   s.t.   Ax + s = b,  s in K           (KKT: Px + q + A'z = 0, s in K, z in K*, s'z = 0)
 * in the form of cuOpt's C API, through its Wolfe dual
 *     min 1/2 x'Px + b'z   s.t.   Px + A'z = -q,  z in K* = R^zero x R+^nonneg x L^d1 x ... x L^dk.
 *
 * cuOpt returns no dual variables for problems with quadratic constraints, so the multipliers z (the
 * contact impulses and bilateral percussions) are made primal variables of the cuOpt problem. The
 * stationarity conditions of the dual reproduce the KKT system above; its optimal value is minus the
 * optimal value of the cone program.
 *
 * cuOpt variables are [x_model, z]. x_model holds the x_j with a nonzero column of P; a variable with
 * a zero column (the percussion of a perfect constraint, whose compliance is zero) must not enter A
 * and has no cost, it is left out and recovered as 0. Every row is an equality (row j of [P A']). The
 * cuOpt objective x'Qx + c'y has no factor 1/2, so Q = P/2. A cone of dimension d >= 2 is the
 * quadratic constraint z_1^2 + ... + z_{d-1}^2 - z_0^2 <= 0 with the head bounded below by 0; a cone
 * of dimension 1 is only that bound.
 */
struct Model {
    CuoptInt num_vars{0};
    CuoptInt num_rows{0};

    // Linear rows in CSR form, all equalities with right-hand side rhs.
    std::vector<CuoptInt> row_offsets;
    std::vector<CuoptInt> columns;
    std::vector<real_t> values;
    std::vector<real_t> rhs;

    // Variable bounds and linear objective.
    std::vector<real_t> lower;
    std::vector<real_t> upper;
    std::vector<real_t> objective;

    // Quadratic objective x'Qx as triplets, each (row, col) once.
    std::vector<CuoptInt> quad_rows;
    std::vector<CuoptInt> quad_cols;
    std::vector<real_t> quad_values;

    // Quadratic cone constraints: cone k has the diagonal entries cone_vars/cone_coeffs in
    // [cone_offsets[k], cone_offsets[k + 1]), head first (coefficient -1), then the tail (+1).
    std::vector<CuoptInt> cone_offsets;
    std::vector<CuoptInt> cone_vars;
    std::vector<real_t> cone_coeffs;

    // Layout: cuOpt variable of x_j (-1 if x_j is left out); z_i is cuOpt variable z_offset + i.
    std::vector<CuoptInt> x_var;
    CuoptInt z_offset{0};

    int numCones() const { return cone_offsets.empty() ? 0 : static_cast<int>(cone_offsets.size()) - 1; }
};

/**
 * @brief Builds the cuOpt model of the cone program (P, q, A, b, dims) into `model`, reusing its
 * buffers. Only the upper triangle of P (including the diagonal) is read; it must be positive
 * semidefinite. Throws std::invalid_argument for inconsistent dimensions, for a variable with a zero
 * column of P that enters A or has a cost (the dual does not determine it), and for sizes beyond the
 * 32-bit indices of cuOpt.
 */
void buildModel(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims, Model& model);

/// x and z of the cone program from cuOpt's primal solution of `model`, and s = b - Ax.
void recoverSolution(const Model& model, const CscMatrix& A, const VectorXr& b, const std::vector<real_t>& primal, VectorXr& x, VectorXr& z, VectorXr& s);

/// Residuals of the KKT conditions of the cone program, each relative to the size of its data.
struct KktResiduals {
    real_t primal{0};           // cone violation of s = b - Ax, / (1 + ||b||_inf)
    real_t dual{0};             // ||Px + q + A'z||_inf / (1 + ||q||_inf)
    real_t dual_cone{0};        // cone violation of z, / (1 + ||z||_inf)
    real_t complementarity{0};  // max over cone blocks |s_j'z_j| / (1 + ||s_j|| ||z_j||)

    real_t max() const;
};

/// KKT residuals of (x, s, z) for the cone program; only the upper triangle of P is read.
KktResiduals kktResiduals(const CscMatrix& P, const VectorXr& q, const CscMatrix& A, const VectorXr& b, const ConeDims& dims, const VectorXr& x, const VectorXr& s, const VectorXr& z);

}  // namespace cardillo::solver::cuopt
