#pragma once

#include "types.hpp"

namespace cardillo::misc {

// How a per-contact 3x3 (normal + 2 tangential) Delassus/Gii block's local rho=(rhoN,rhoT,rhoT)
// step-size/preconditioner scale is derived. Mirrors Siconos's own three strategies
// (fc3d_AlartCurnier_functions.c: compute_rho_split_spectral_norm[_cond]/compute_rho_spectral_norm)
// -- checked directly against that codebase, not guessed. Shared between the local semismooth
// Newton solve (solver/local_contact_newton.cpp) and CondensedAssembler's GiiInv construction
// (assembly/condensed_assembler.cpp), both of which need the same "how strongly is this block's
// own normal/tangential response coupled" answer for the same 3-row contact block:
//   Split: rhoN=1/G(0,0), rhoT=1/lambda_max(G_TT) (G_TT the 2x2 tangential sub-block) -- treats
//     normal and tangential as decoupled 1- and 2-dof sub-problems for the step-size choice only.
//     Matches Siconos's compute_rho_split_spectral_norm.
//   FullSpectral: a single rho = 1/lambda_max(sym(G)) applied to all three components, from the
//     FULL 3x3 block's largest eigenvalue (symmetrized -- G need not be exactly symmetric once an
//     implicit-gyroscopic body is involved). Matches Siconos's compute_rho_spectral_norm.
enum class ContactRhoStrategy { Split, FullSpectral };

// Returns (rhoN,rhoT,rhoT). Any component is -1 if the corresponding sub-block is (near-)singular
// (<= eps) -- callers must check for a non-positive component and treat that as "this G is not
// usable for a rho-based step size," falling back to whatever they'd otherwise do.
Vector3r computeContactRho(const Matrix33r& G, ContactRhoStrategy strategy, real_t eps = (real_t)1e-12);

}  // namespace cardillo::misc
