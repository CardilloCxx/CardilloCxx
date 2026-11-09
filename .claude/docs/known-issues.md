# Known open issues

Status as of 2026-10-31 (after the build restructure and the bug-fix PR of the cleanup run). Check
whether an item was fixed before re-investigating; remove items once they are resolved.

## Scenes

- **hangbride** (`hangbridge` target): diverges with the condensed solver under the position-level
  law; θ = 0.5 in its config cannot work for perfect constraints, θ = 0.7 finished once but with a
  suspicious kinetic energy (~8.6e4). With PGS it diverges near step 245 independently of recent fixes.
- **ciliaWalker** (`cilia_walker`): the bidirectional and ring walkers (thin rectangular cilia, 12
  elements, h = 1e-4) diverge in the condensed Gauss-Seidel sweep after ~0.1 s unless
  `condensed.true_schur=true` (or h = 5e-5). The Coal crash that followed is fixed by the monotone-chain
  convex hull in the Coal fork (see `tests/test_convex_hull.cpp`).
- **cable_bunch**: the segfault inside Coal `computeSupportSetConvexHull` (box-box patch, around step
  1247) came from the fork's convex hull (non-strict sort comparator, out-of-bounds reads on duplicate
  points); fixed in the fork, not yet rerun to step 1247 (several hours at 4000–10000 iterations/step).
- **chainFountain** (Clarabel, dt = 2e-5): Clarabel fails around step 2428.
- **conveyor_merge**: waterfall height 0 (and 0.1 mm) gets the parcel stuck on the coplanar box seam of
  the merge belt (horizontal contact normal); works from 1 mm. At 200 Hz the 10 mm case is not
  converged in yaw (−5.2° vs −3.9° at 1 kHz).
- Scenes with random initial data (`billiard`, `parcel`, `constraintTest`: `random_device`/`UnitRandom`)
  are not reproducible run to run.

## Code

- Velocity setters still used in scenes: `parcel` (belt as 1e8 kg body + `setLinearVelocity`), 3DPrinter
  (`setConstraintAngularVelocity`). Goal: pose-only drivers / surface velocity, then remove the setters.
- `C_BeamHull`/`C_Collider_BeamHull` store radii as `float` (with `(float)` casts in `BeamFactory`)
  while the engine uses `real_t`.
- PGS `DinvDiag` inverts the contact rows per row, i.e. scales the tangential directions
  anisotropically; its fixed point satisfies a distorted Coulomb law (condensed uses a common
  tangential step length).
- PGS ignores gyroscopic forces with `moreau.implicit_gyroscopy=true`; condensed ignores restitution.
- Clarabel exports zero tangential percussions in the VTK output (not investigated; test solver).
- `docs/chapters/using_the_engine/beams.rst` ("Rod model") cites the gitignored paper.
- QOCO fails ("Numerical Error") on the rigid pendulum of `test_position_level` (perfect constraint, no
  contacts), so that test fails for QOCO.
- MOSEK: the pendulum percussion of `test_position_level` deviates by ~4e-6 (relative) from projected
  Jacobi independent of the tolerances (suspect: epigraph reformulation of ½xᵀPx); MOSEK is therefore
  not in that test.
- Default scene configs that select an optional conic backend (`chainFountain`, `jenga`, `painleve`,
  `strandbeest`, `tennis`) only run when that backend is built.
- `docs/chapters/moreau_time_stepping.rst` ("Scaled Block Linear System") still uses the old sign
  convention Lambda = -dt lambda; the code and theory manual use pi = h lambda acting as +W pi.
