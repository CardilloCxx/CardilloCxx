#pragma once

#include "../../misc/types.hpp"

namespace cardillo::solver {

struct NewtonACParams {
    int maxIters{8};
    real_t tol{(real_t)1e-10};
};

// Guarded Alart-Curnier semismooth Newton solve for one frictional (3-row) contact block, in the
// same internal sign convention as ProjectedGaussSeidel/CondensedSolver: lambda[0] (normal) <= 0,
// Coulomb disk radius = -mu*lambda[0]. See local_contact_newton.cpp for the derivation.
//
// G: the block's 3x3 mass-only Delassus matrix (RowBlock::Gii; contacts carry zero compliance, so
//    this is the full local matrix -- no separate "+compliance" step needed).
// r: the block's current residual, r = rhs_row - Ja*u_corr_A - Jb*u_corr_B, evaluated at the
//    CURRENT `lambda` (i.e. blockResidual()'s return value for this block).
// lambda: in/out. On success, overwritten with the converged impulse (lambda[0] <= 0 guaranteed by
//    construction of every case branch). Left UNCHANGED on failure.
// Returns false if the generalized Jacobian is (near-)singular, a Newton step is non-finite, a
// step fails to reduce ||Phi|| after the first iteration, or maxIters is exhausted without
// meeting tol -- caller should fall back to a plain projection step in all of these cases.
bool solveContactBlockNewtonAC(const Matrix33r& G, const Vector3r& r, real_t mu, Vector3r& lambda, const NewtonACParams& params);

}  // namespace cardillo::solver
