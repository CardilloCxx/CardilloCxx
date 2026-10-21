#pragma once

#include "../../misc/contact_rho.hpp"
#include "../../misc/types.hpp"

namespace cardillo::solver {

// See misc/contact_rho.hpp for the Split/FullSpectral strategies and their theory -- shared with
// CondensedAssembler's GiiInv construction, which needs the same per-block rho estimate.
using NewtonRhoStrategy = misc::ContactRhoStrategy;

struct NewtonACParams {
    int maxIters{8};
    real_t tol{(real_t)1e-10};
    NewtonRhoStrategy rhoStrategy{NewtonRhoStrategy::Split};
};

// Guarded Alart-Curnier semismooth Newton solve for one frictional (3-row) contact block (theory
// manual, sec:local_newton), in the negative variant lambda = -proj_C(rho*xi - lambda) of all
// projection-based solvers: lambda[0] (normal) >= 0, Coulomb disk radius = mu*lambda[0] (see
// misc/contact_projection.hpp). See local_contact_newton.cpp for the derivation.
//
// G: the block's 3x3 mass-only Delassus matrix (RowBlock::Gii; contacts carry zero compliance, so
//    this is the full local matrix -- no separate "+compliance" step needed).
// xi: the block's current row velocity, xi = xi^free + Ja*du_A + Jb*du_B, evaluated at the
//    CURRENT `lambda` (i.e. blockResidual()'s return value for this block).
// lambda: in/out. On success, overwritten with the converged impulse (lambda[0] >= 0 guaranteed by
//    construction of every case branch). Left UNCHANGED on failure.
// Returns false if the generalized Jacobian is (near-)singular, a Newton step is non-finite, a
// step fails to reduce ||Phi|| after the first iteration, or maxIters is exhausted without
// meeting tol -- caller should fall back to a plain projection step in all of these cases.
bool solveContactBlockNewtonAC(const Matrix33r& G, const Vector3r& xi, real_t mu, Vector3r& lambda, const NewtonACParams& params);

}  // namespace cardillo::solver
