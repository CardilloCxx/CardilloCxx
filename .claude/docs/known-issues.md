# Known open issues

Status as of 2026-10-20 (after the moreau_theta_position branch was merged as PRs 1–13). Check
whether an item was fixed before re-investigating; remove items once they are resolved.

## Scenes

- **hangbride** (`hangbridge` target): diverges with the condensed solver under the position-level
  law; θ = 0.5 in its config cannot work for perfect constraints, θ = 0.7 finished once but with a
  suspicious kinetic energy (~8.6e4). With PGS it diverges near step 245 independently of recent fixes.
- **3DPrinter**: since #77 (node-based `InterSegmentHull` sampling) the belt needs ~800 condensed
  iterations per step instead of ~150. About half comes from shear prestress at the pulleys (node
  frames follow the spline tangent while `relaxedInStraightConfiguration()` zeroes the shear rest
  strain); the rest from the node layout itself. Not caused by the pretension. Topic for the #77 author.
- **gears**: aborts with "PjAssembler::solveS called but S matrix is not factorized" (also on main
  before the position-law migration).
- **cable_bunch**: segfault inside Coal `computeSupportSetConvexHull` (box-box patch) around step 1247.
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
- Default `condensed.projection_rho_strategy=diagonal` scales the tangential directions anisotropically;
  its fixed point satisfies a distorted Coulomb law (`split`/`full` are consistent).
- PGS ignores gyroscopic forces with `moreau.implicit_gyroscopy=true`; condensed ignores restitution.
- Clarabel exports zero tangential percussions in the VTK output (not investigated; test solver).
- `docs/chapters/using_the_engine/beams.rst` ("Rod model") cites the gitignored paper.
