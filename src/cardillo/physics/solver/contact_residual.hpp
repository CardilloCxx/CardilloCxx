#pragma once

#include "../../misc/types.hpp"
#include "../assembly/dynamics_assembler.hpp"

namespace cardillo::solver {

// Solver-independent accuracy check of the contact part of a Moreau-θ step: the blockwise
// natural-map residual of the stored contact percussions `Contact::last_impulse` against the
// contact velocities xi = W u + offset of the new velocities `u`, with contact-local step sizes
// from the Delassus blocks (split strategy). Evaluated once, nothing is updated. A value ≤ 1 means
// that every contact satisfies the tolerance (tol_abs, tol_rel). Returns NaN for non-finite data.
// Must be called after the solve and before the new velocities are written to the world.
real_t contactNaturalMapResidual(physics::DynamicsAssembler& dyn, const VectorXr& u, real_t tolAbs, real_t tolRel);

}  // namespace cardillo::solver
