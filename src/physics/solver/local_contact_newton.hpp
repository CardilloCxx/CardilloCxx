#pragma once

#include "../../misc/types.hpp"

namespace cardillo::solver {

// How the per-contact Newton step's diagonal preconditioner rho=(rhoN,rhoT,rhoT) is derived from
// the block's own local Delassus matrix G. Mirrors Siconos's own three strategies
// (fc3d_AlartCurnier_functions.c: compute_rho_split_spectral_norm[_cond]/compute_rho_spectral_norm)
// -- checked directly against that codebase, not guessed:
//   Split (default, this codebase's original/only strategy): rhoN=1/G(0,0), rhoT=1/lambda_max(G_TT)
//     (G_TT the 2x2 tangential sub-block) -- treats normal and tangential as decoupled 1- and
//     2-dof sub-problems for the STEP-SIZE choice only (the full 3x3 G is still used in the
//     Newton residual/Jacobian either way). Matches Siconos's compute_rho_split_spectral_norm.
//   FullSpectral: a single rho = 1/lambda_max(sym(G)) applied to all three components, from the
//     FULL 3x3 block's largest eigenvalue (symmetrized -- G need not be exactly symmetric once an
//     implicit-gyroscopic body is involved) -- captures normal-tangential coupling in the
//     step-size choice, which Split ignores. Matches Siconos's compute_rho_spectral_norm.
enum class NewtonRhoStrategy { Split, FullSpectral };

struct NewtonACParams {
    int maxIters{8};
    real_t tol{(real_t)1e-10};
    NewtonRhoStrategy rhoStrategy{NewtonRhoStrategy::Split};
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
