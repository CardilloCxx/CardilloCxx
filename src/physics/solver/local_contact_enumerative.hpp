#pragma once

#include "../../misc/types.hpp"

namespace cardillo::solver {

struct EnumerativeParams {
    // Generic near-zero / degeneracy guard, reused for every threshold this solver needs: |wN| too
    // small to form the sliding case's Schur complement, |Delta| too small to trust the sliding
    // reconstruction, alpha required strictly > eps (not just >= 0), and the disk/cone membership
    // checks' slack for take-off/sticking validity. Deliberately a single scalar -- this solver has
    // no iteration count or convergence tolerance to tune (it is a closed-form case enumeration,
    // not an iterative method), so there is nothing else to expose here.
    real_t eps{(real_t)1e-10};
};

// Enumerative (Bonnefon & Daviet, "Quartic formulation of Coulomb 3D frictional contact", INRIA
// RT-0400 2011; used as in Daviet, Bertails-Descoubes & Boissieux, "A Hybrid Iterative Solver for
// Robustly Capturing Coulomb Friction in Hair Dynamics", SIGGRAPH Asia 2011) guarded local solve
// for one frictional (3-row) contact block, in the same internal sign convention as
// ProjectedGaussSeidel/CondensedSolver/solveContactBlockNewtonAC: lambda[0] (normal) <= 0, Coulomb
// disk radius = -mu*lambda[0].
//
// Unlike solveContactBlockNewtonAC, this is not an iterative correction from the current lambda --
// it enumerates the three disjunctive Coulomb cases (take-off / sticking / sliding) directly, the
// sliding case reducing to a scalar quartic (see quartic_solver.hpp and local_contact_enumerative.cpp
// for the derivation and the exact sign/notation translation from the paper's r_N>=0 convention).
// It still folds in the incoming `lambda` (via the same q = r + G*lambda construction
// solveContactBlockNewtonAC uses at local_contact_newton.cpp:118), since `r` in this codebase is a
// residual relative to the current lambda, not an absolute free term -- see this function's first
// line.
//
// G, r: exactly as documented in local_contact_newton.hpp (blk.Gii, blk.dim==3 block residual). G
//    need not be symmetric (moreau.implicit_gyroscopy) -- see the .cpp file's derivation, which is
//    verified correct for this case, not just the symmetric one the source paper assumes.
// lambda: in/out, same contract as solveContactBlockNewtonAC -- read for the q=r+G*lambda fold-in,
//    overwritten with the accepted solution on success (lambda[0] <= 0 guaranteed by construction),
//    left UNCHANGED on failure.
// Returns false if no case (take-off/sticking/sliding) yields a valid candidate. For a symmetric
// positive-definite G this should essentially never happen (the three disjunctive cases jointly
// cover the whole NCP); for a non-symmetric G (the gyroscopic case) it is a real, expected
// possibility -- existence/uniqueness of the Coulomb NCP's solution is only guaranteed for
// symmetric positive-definite local Delassus matrices. Symmetrically, when MULTIPLE candidates
// validate at once (also only possible for non-symmetric G), the one minimizing a simple
// Alart-Curnier fixed-point residual is returned -- see candidateResidual() in the .cpp file.
bool solveContactBlockEnumerative(const Matrix33r& G, const Vector3r& r, real_t mu, Vector3r& lambda, const EnumerativeParams& params);

}  // namespace cardillo::solver
