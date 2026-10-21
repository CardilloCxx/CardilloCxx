# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

CardilloCxx is a C++20 nonsmooth multibody engine: rigid bodies, point masses, Cosserat rods (beams) and
soft bodies with frictional contact, integrated by a Moreau-θ time-stepping scheme with several
solver backends. Topic notes with more detail live in `.claude/docs/`:

- `.claude/docs/time-integration.md` – position-level force law, sign conventions, drivers, joints
- `.claude/docs/collision.md` – Coal fork, contact patches, beam hulls, mesh frames, surface velocity
- `.claude/docs/known-issues.md` – open problems that are known and need not be re-investigated

## Commands

```bash
cmake -S . -B build                       # default build type RelWithDebInfo, -march=native, tests ON
cmake --build build -j$(nproc)            # library, all scene executables, tests
cmake --build build --target chain        # a single scene (target names: examples/CMakeLists.txt)
cd build && ctest                         # all tests; single test: ctest -R test_torsional_oscillator
build/tests/test_position_level           # tests are plain executables, run them directly for output
cd docs && uv run make html               # Sphinx + Doxygen docs
```

- Run scenes **from the repository root** (mesh paths like `res/meshes/...` are relative to the cwd):
  `build/examples/<target> examples/scenes/<dir>/scene.config`. (The README's `build/bin/` is outdated.)
- Short runs: copy the config and append overrides; a key given later wins, e.g.
  `printf "\nsim.T=0.05\noutput.interval_steps=0\n" >> copy.config`. `.claude/scripts/run_scenes.sh`
  automates this for several scenes (see the `scene-smoke` skill). Full scenes can take minutes to hours.
- Coal is fetched by FetchContent from a fork (see `.claude/docs/collision.md`). In extra worktrees,
  reuse an existing checkout with `-DFETCHCONTENT_SOURCE_DIR_COAL=<path to coal-src>`.
- Clarabel needs Rust/cargo. Its Rust core is built in cargo's release profile for Release,
  RelWithDebInfo and MinSizeRel (`src/CMakeLists.txt` presents them as Release to Clarabel's CMake);
  only Debug gives a debug Rust build (15–60× slower solves). Check with
  `ldd build/src/libcardillo.so | grep clarabel` (path must contain `target/release`) before timing.
- Changing `config::Config` (or other headers used by `example_main.hpp`) requires rebuilding **all**
  scene executables; stale binaries crash with a Config layout mismatch that looks unrelated.
- Formatting: `.clang-format` (Google based, 4 spaces, column limit 200). `.clang-tidy` checks unused
  code only. A clean build has no compiler warnings; keep it that way.

## Architecture

**Facade and ownership.** `physics::PhysicsEngine` (`src/physics/api`) is the scene-facing API. It
owns the `World` (EnTT registry with ECS components from `src/physics/ecs_types.hpp`, mesh assets,
constraint patterns), `collision::CollisionCoal` and the `PhysicsPipeline`, which owns the
`DynamicsAssembler`, the solver selected by `solver.type`, and the `MoreauIntegrator`.

**One step** (`src/physics/integration/moreau.cpp`): drift positions by (1−θ)h → advance prescribed
constraint rates → `DynamicsAssembler::updateStateDependentTerms` (kinematic drivers evaluated at
t_{n+θ}, rigid-state cache, beam element geometry, Coal collision → contacts, contact Jacobian `W`,
constraint rows `Wg`/`Wγ`, compliances, position errors) → `solver.solve(dt, θ)` returns the new
velocities → drift by θh.

**DynamicsAssembler** (`src/physics/assembly`) is the single source of the discretization: flat
velocity vector with per-body offsets (`m_body_vel_offsets`; rigid bodies 6 dofs, point masses 3),
contact rows (each contact side may have several weighted attachments, used by beam hulls),
spring/damper rows from all `ConstraintPattern`s, and the Moreau-θ terms `springBias()`,
`damperBias()`, `springComplianceScale()`, `damperComplianceScale()`. Every solver must use these
instead of its own copies. Each solver has a matching assembler (`pj_`, `condensed_`, `pgs_`,
`qoco_`, `clarabel_`, `conicxx_assembler`).

**Solvers** (`src/physics/solver`). The production solvers are **projected Jacobi** (`pj`: direct
factorization of the bilateral saddle-point system, projected iteration on contacts) and
**condensed** (block-sparse Schur complement; `condensed.true_schur=true` eliminates the bilateral
subsystem exactly, needed for stiff rods). Test and validate mainly with these two. PGS, CG,
QOCO, Clarabel and ConicXX are convexified/test backends.

**Constraints** (`src/physics/constraints`): `getConstraint()` returns a `ConstraintResult`
(Jacobian columns `WgA/WgB`, compliances `Crows/Arows`, `positionError`, prescribed rate). Main
patterns: `LinearDistanceConstraint`, `TranslationRotationConstraint` (all joints; hinge = x-axis),
`BeamConstraint` (Cosserat rod element between two nodes, Petrov–Galerkin force directions).

**Bodies** (`src/physics/assets`): `RigidBodyFactory`, `BeamFactory` (rods = chains of rigid nodes;
colliders `RigidBodyPrimitive` or `InterSegmentHull`), soft bodies, point masses. Dynamic mesh bodies
live in the principal-axes frame of their mesh (see collision notes).

**Config** (`src/config`): `key = value` files, `#` comments. Unknown keys warn with `file:line`.
Adding a key means: field with default in `config.hpp`, parser branch in `config.cpp`, row in the
table of `docs/chapters/using_the_engine/setup.rst`. Scene-specific values use `scene.<name> = value`
and are read via `engine.config().sceneParam("<name>", fallback)`.

**Scenes** (`examples/scenes/<dir>/`): a `SceneBase` subclass with `populate(engine)` and optionally
`updateScene(engine, t)` (step size: `engine.timeStep()`; the old overload with `dt` is deprecated),
registered with `cardillo_add_scene_example(...)` in `examples/CMakeLists.txt`; `examples/main.cpp` +
`example_main.hpp` run it with a config file. Output: VTK (`output.folder`, `output.interval_steps`;
0 disables) and a CSV of bodies tagged with `engine.track(e, name)`.

**Tests** (`tests/`, ctest): `test_joint_jacobian`, `test_beam_jacobian` (W vs. d/dt g),
`test_beam_rest_configuration`, `test_position_level` (all solvers), `test_torsional_oscillator`
(closed-form solution, convergence order), `test_condensed_true_schur` (condensed vs. PJ).

## Design decisions (do not undo without discussion)

- **Position-level Moreau-θ force law only**: no Baumgarte stabilization, no `lambda_theta` variant, no
  velocity-level law. Perfect constraints need θ > 1/2 (warning otherwise). Details:
  `.claude/docs/time-integration.md`; theory in `docs/chapters/moreau_time_stepping.rst` and, if present
  locally, the paper `LaTeX/theory_doc.tex` (gitignored, separate repository).
- **Kinematic drivers are pose-only**: `addTrajectory(e, poseFunc)`; the velocity is derived (mean over
  the step). Prescribe velocities through their integrated pose. Remaining `setLinearVelocity`/
  `setAngularVelocity`/`setVelocityByForce` uses in scenes are to be removed, not extended.
- **Compliances, not stiffnesses**, in `Crows/Arows`: stiffness ∞ → compliance 0 (perfect constraint),
  stiffness ≤ 0 → compliance ∞ (free row, not assembled). No magic thresholds.
- **Sign conventions** (theory manual, `LaTeX/theory_doc.tex`): `Contact::last_impulse = (pn, pt1, pt2)`
  acts on body B (normal from A to B), `pn ≥ 0`, in every solver. Bilateral percussions are
  π = h·λ (`Lambda_g`, `Lambda_gamma` hold π_g, π_γ). Every block acts as +W·π with the force
  directions W of `DynamicsAssembler` unchanged; u = u_free + M⁻¹Wπ, ξ = Wᵀu + Ĉπ + b. Contacts use the
  negative-variant projection λ = −proj_C(ρξ − λ) (`misc/contact_projection.hpp`, shared with PJ).
- **No conversions to Coal**: `real_t == coal::Scalar`, `Vector3r == coal::Vec3s` (checked once in
  `src/collision/coal_types.hpp`); our own code uses `Vector3r`, never `coal::Vec3s`, and no casts.
- Torsional springs belong on the joint x-axis (exact unwrapped twist angle; y/z rows measure sines).

## Code style expected here

- No `(void)x`: leave unused parameters unnamed (`real_t /*dt*/`), delete unused locals.
- No hard-coded "magic" numbers or hacks to make a scene run; find the actual cause.
- Every scalar/vector/matrix member of an ECS component has a default member initializer.
- Prefer fixed-size Eigen types; compute shared quantities once and pass them (e.g. joint triads), copy
  small vectors into locals in hot loops (measured faster than repeated column blocks).
- Comments describe the current code, not its history or removed variants.

## Workflow

- Do not commit or push unless asked; confirm before force-pushes. Larger work is split into small,
  reviewable PRs (stacked branches `pr/NN-topic`), each building warning-free with all tests passing.
- PR descriptions: title separately, body as GitHub markdown inside a ```markdown block.
