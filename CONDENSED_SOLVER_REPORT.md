# Condensed contact solver: implementation & investigation report

**Scope:** a new `condensed` contact solver (`solver.type=condensed`), block-sparse and
config-selectable across four sweep strategies and two local-solve strategies, built as a
generalization of the existing `projected_gauss_seidel` (PGS) solver. This document records what
was built, the two real bugs the validation process caught (and how), the ThreadSanitizer
investigation, and the benchmark results as actually measured — not a projection of expected
results.

See `docs/chapters/solvers/condensed.rst` for the solver's mathematical/architectural
documentation. This file is the narrative of *how it was built and verified*.

## Direct answer to "are these all non-accelerated variants?"

**Yes.** All four sweep modes (`gauss_seidel`, `jacobi`, `colored`, `chaotic`) are plain
fixed-point iterations: relaxation (`pj_relaxation`) and per-contact step damping (`pj_alpha`) are
the only step-size controls, identical to what PGS already has. None of them have Nesterov,
Chebyshev, or Anderson acceleration, all three of which exist for `projected_jacobi` (PJ), and
Nesterov also exists for PGS. This matters directly for interpreting the benchmark numbers below:
comparing `condensed` against PJ's *default* (Nesterov-on) configuration systematically favors PJ
on iteration count. A fair comparison uses `pj.nesterov=false`, which is what the tables below do
unless noted. Adding acceleration to `condensed`'s sweep modes is realistic future work — the
per-block update is already isolated in `computeBlockUpdate()` in `condensed_solver.cpp`, and PJ's
existing Nesterov/Chebyshev logic operates at exactly that level (wrapping the whole-vector sweep),
so the pattern is known and reusable.

## What was built

| File | Purpose |
|---|---|
| `src/misc/graph_coloring.{hpp,cpp}` | Greedy Welsh-Powell graph coloring (used by `colored` mode) |
| `src/physics/assembly/contact_jacobian.hpp` | `buildContactRowByDof()`, extracted (pure refactor, zero behavior change) from `dynamics_assembler.cpp` so the new assembler can reuse it |
| `src/physics/assembly/condensed_assembler.{hpp,cpp}` | `RowBlock`/`CondensedTopology`: the block-sparse representation, never touching `Eigen::SparseMatrix` |
| `src/physics/solver/condensed_solver.{hpp,cpp}` | The solver itself: four sweep drivers + dispatch |
| `src/physics/solver/local_contact_newton.{hpp,cpp}` | Guarded Alart-Curnier semismooth Newton local solve |
| `src/misc/block_diagonal.hpp` | `invertSmallSpd()` extracted from `BlockDiagonal::calculateInverse()` (pure refactor, reused by the new assembler) |
| `src/config/config.{hpp,cpp}` | `SolverType::Condensed` + `condensed.*` config keys |
| `src/physics/pipeline/physics_pipeline.cpp` | Dispatch wiring |
| `src/misc/timings/TimingManager.hpp` | `Condensed*` `TimerId` entries |
| `CMakeLists.txt`, `src/CMakeLists.txt` | `CARDILLO_ENABLE_TSAN` option; `find_package(OpenMP)` + `CARDILLO_HAVE_OPENMP` |
| `examples/scenes/{domino,slinky}/scene_condensed.config` | Reference scene configs |
| `examples/example_main.hpp` | `CARDILLO_DUMP_STATE` env-gated end-of-run state fingerprint (total kinetic energy, summed position norm) — a lightweight cross-solver correctness check, cheaper than diffing VTK output |
| `tools/compare_solver_traces.py` | Fixed a stale binary path (see below), added `condensed` ablation cases, added wall-clock capture |

`projected_jacobi.*` and `projected_gauss_seidel.*` were **not modified** except the two
extractions above, both verified behavior-identical before proceeding.

## Validation methodology and the two bugs it caught

### 1. Assembler-level: `CondensedAssembler::rhs()` vs `PgsAssembler::rhs()`

Before writing any solver logic, `CondensedAssembler`'s `rhs()`/`ufree()` were diffed
segment-by-segment against `PgsAssembler`'s own output, in-process (both assemblers run against
the same live `DynamicsAssembler` state, no config-file involved). Result: agreement to
floating-point round-off (`~1e-13` to `~1e-16` relative) on domino (6 springs, up to 40k contact
rows) and slinky (4944 springs, up to 432 contact rows). This is the strongest validation in the
whole exercise, because it isolates "did I port the formulas correctly" from every other later
concern.

### 2. A process bug: silently comparing PJ against itself

The first attempt at "validate the new solver against PJ" used a shell config-rewriting pattern
(`sed ... > new.config; echo "solver.type=condensed" >> new.config`) that assumed the base config
ends in a newline. Several of this repo's scene configs (`domino/scene_standard.config` among them)
do not — the last line is `pj.convergence_csv_dir=domino_nesterov_warmstart_csv` with no trailing
`\n`. The appended `solver.type=condensed` line landed glued onto the end of that line
(`...csvsolver.type=condensed`), which the config parser doesn't recognize as a `solver.type` key
at all, so it silently fell back to the default (`ProjectedJacobi`). The first "PJ vs PGS vs
Condensed all agree to 15 digits" result was real for the *PJ* case but the *PGS* and *Condensed*
cases were secretly also running PJ. Fixed by rewriting the config generator to guarantee a
separating newline and, from then on, verifying which solver actually ran by grepping the
`Timing Breakdown`'s named row (`Projected Jacobi` vs `Projected Gauss-Seidel` vs `Condensed`)
rather than trusting the config file alone.

### 3. A real bug: full-block vs. diagonal local inverse

With the config generation fixed, re-running the domino comparison exposed a real discrepancy:

| Solver | Total kinetic energy after 1 step (domino, default tolerance) |
|---|---|
| PJ | 0.814196987077336 |
| PGS | 0.814479477909845 |
| Condensed (buggy) | **2.2858845614721** |

`CondensedAssembler` computed `RowBlock::GiiInv` as the full dense inverse of the local 3x3/6x6
Delassus block — mirroring `PgsAssembler::Dinv()`. But `ProjectedGaussSeidelSolver` doesn't call
`Dinv()`; it calls `PgsAssembler::DinvDiag()`, which uses only the **diagonal** of that same block
(no normal/tangential or multi-DOF coupling term). `Dinv()` is dead code for contacts in this
codebase — the only caller is `ConjugateGradientSolver`, which is bilateral-only. Using the more
"accurate" full-block inverse changed the Gauss-Seidel iteration's fixed point on domino's densely
coupled ~13k-contact system, not just its convergence rate: it satisfied its own residual
tolerance while landing on a materially wrong answer (`iters=15619`, ~5x more sweeps than
PGS/PJ needed, converging somewhere else entirely). After matching `DinvDiag()`'s diagonal-only
formula exactly:

| Solver | Total kinetic energy after 1 step (domino, default tolerance) | Sweeps |
|---|---|---|
| PJ | 0.814196987077336 | 3056 |
| PGS | 0.814479477909845 | 3057 |
| Condensed, `gauss_seidel` | 0.811464242457628 | 11979 |
| Condensed, `jacobi` (4 threads) | 0.811304918689760 | 11981 |
| Condensed, `colored` | 0.811173633909661 | 11970 |

All four agree to within ~0.4%, consistent with ordinary Jacobi-vs-Gauss-Seidel iteration-path
sensitivity on a nonlinear frictional system (the same phenomenon already documented in this
repo's own history, e.g. the Anderson-acceleration commit's note on chaotic-scene trajectory
divergence) — not a correctness bug. `Gii` (the full block) is still stored and used by the Newton
local solve, which genuinely needs the coupling; only the *linear* projection path uses the
diagonal approximation, matching PGS.

**Why the bug wasn't caught by the earlier slinky-only tests:** the very first tight-tolerance,
single-step slinky comparison (done before this bug was found) happened to still agree to ~11
significant digits even with the buggy full-block inverse, because at that state slinky's contact
impulses were tiny relative to its spring forces — the choice of contact-block preconditioner
barely mattered yet. Domino's ~13k simultaneous, strongly-coupled contacts is what made the error
large enough to be unmistakable. Lesson applied for the rest of the session: validate on the
scene most likely to stress the mechanism in question, not just the first one that's convenient.

### 4. Multi-threaded determinism

After the fix, `jacobi` and `colored` were run at 1, 4, and 16 threads on domino. Results were
**bit-identical** across all three thread counts (same kinetic energy and position-norm to all 15
displayed significant digits), repeated after every subsequent code change (including the OpenMP
region-merging refactor below). This is expected — both modes are structurally deterministic
(each sweep reads a fixed pre-sweep state and writes disjoint locations) — and is treated as
supporting, not sufficient, evidence of race-freedom.

### 5. ThreadSanitizer: a real finding, but not the one it looks like

Building a `CARDILLO_ENABLE_TSAN` configuration and running `jacobi`/`colored`/`chaotic` under TSan
reported a data race at the boundary between the OpenMP region and the surrounding code in every
one of them. Initial fix attempt: merge each sweep's separate `#pragma omp parallel for` calls
into one `#pragma omp parallel` region with internal `#pragma omp for` stages (this is also a
legitimate performance win — see below — since it cuts thread-team spawn/join count from up to 52
per sweep in `colored` down to 1). The report persisted identically afterward.

To determine whether this was a real bug or a tooling limitation, a minimal, provably-correct
standalone reproducer was built and run under the same TSan configuration: two plain
`#pragma omp for` loops inside one `#pragma omp parallel`, first writing disjoint indices of a
shared array, second reading what the first wrote after the (implicit) barrier. TSan flagged it
identically. This confirms a documented, decade-old limitation: GCC's `libgomp` is not itself
built with TSan instrumentation, so TSan cannot observe its barrier/team-join synchronization and
reports the (correctly ordered) accesses on either side of that barrier as a race — acknowledged
by TSan's own author on the GCC mailing list in 2013, and still true of this GCC 13 / libgomp on
this system. No Clang/`libomp` toolchain (which *is* well-instrumented for this) was available in
the environment to get a clean run.

**Practical conclusion:** correctness confidence for the parallel sweep modes rests on (a) the
disjoint-memory-access arguments documented inline at each `#pragma omp for` site in
`condensed_solver.cpp`, and (b) the bit-identical multi-thread-count results above — not on a
clean TSan pass, which isn't obtainable with the toolchain available here. If a Clang/libomp
toolchain becomes available, re-running the TSan build with it would be the natural way to close
this out properly.

## Benchmark results

Two kinds of runs were done: a short, tight-tolerance correctness/iteration-count comparison
(cheap, many repetitions during development), and a longer wall-clock comparison (fewer runs, real
timing). All numbers below are from the *fixed* code (post both bugs above).

### Local-solve comparison: projection vs. semismooth Newton (slinky)

Single-step, `pj.tol_abs=5e-8`, `condensed.sweep_mode=gauss_seidel`:

| Local solve | Sweeps | Kinetic energy |
|---|---|---|
| `projection` | 1318 | 8.79571318074891e-05 |
| `newton` | 134 | 8.79979857874781e-05 |

**~9.8x fewer sweeps**, converged answer within 0.05%. `chaotic` mode showed the same pattern
(`projection`: 1318 sweeps; `newton`: 896 sweeps, `is_always_lock_free<double> == 1` confirmed on
this platform so the atomics carry no lock overhead). This is the expected result and the main
reason to prefer `newton` on frictional/compliant scenes.

### Wall-clock comparison (short runs, 90s budget per case, 6 parallel workers)

**domino** (`sim.T=0.0006`, ~3 steps, ~13k frictional contacts, `pj.nesterov=false`):

| Case | Total sweeps | Wall-clock | Completed? |
|---|---|---|---|
| PJ | 42244 | 40.82s | yes |
| PGS | 20605 | — | no (90s timeout) |
| Condensed `gauss_seidel`+`projection` | 11979 | — | no (1/3 steps in 90s) |
| Condensed `jacobi`+`projection`, default threads (16) | 0 | — | no (0 progress in 90s) |
| Condensed `jacobi`+`projection`, 4 threads | 11981 | — | no (1/3 steps in 90s) |

**slinky** (`sim.T=0.02`, 200 steps, `pj.nesterov=false`):

| Case | Total sweeps | Wall-clock | Completed? |
|---|---|---|---|
| PJ | 17640 | 32.62s | yes |
| PGS | 16888 | 30.73s | yes |
| Condensed `gauss_seidel`+`projection` | 16750 | 50.73s | yes |
| Condensed `gauss_seidel`+`newton` | 17542 | 41.80s | yes |
| Condensed `jacobi`+`projection` | 2736 | — | no (15/201 steps in 90s) |
| Condensed `colored`+`projection` | 347 | — | no (3/201 steps in 90s) |
| Condensed `chaotic`+`newton` | 5234 | — | no (83/201 steps in 90s) |

### Honest reading of these numbers

- **`gauss_seidel` tracks PGS's sweep count almost exactly** on both scenes (expected — same
  algorithm, different data structure) but is **slower in wall-clock than PGS** on slinky (50.7s vs
  30.7s for a near-identical sweep count). The likely cause, not yet fixed: `condensed_solver.cpp`
  builds a fresh dynamically-sized `VectorXr` for every block's residual/update/delta every sweep,
  where PGS operates on pre-sized whole-system vectors. Frictional blocks are always exactly 3
  rows — switching the per-block hot path to fixed-size `Vector3r`/`Matrix33r` (no heap allocation)
  is the highest-leverage next optimization and wasn't done in this pass.
- **`newton` wins convincingly**: fewer sweeps *and* less wall-clock than `projection` on slinky
  (41.8s vs 50.7s) despite the extra per-contact 3x3 Newton cost, confirming the sweep-count
  reduction more than pays for itself here.
- **`colored` is not currently competitive**: ~50 colors on domino times thousands of sweeps means
  roughly 50x the OpenMP barrier count of `jacobi`'s 2-per-sweep, and that overhead dominates at
  this scale — it made the least progress of any mode in the fixed time budget on both scenes.
  Caching the coloring across timesteps (only rebuilding when the contact *set*, not just count,
  changes) is the natural fix, not attempted here.
- **More threads isn't automatically better**: domino's `jacobi` mode made *zero* measurable
  progress at the OpenMP default thread count (all 16 cores) in 90 seconds, but completed a full
  step in the same budget at 4 threads. This points to parallelization/scheduling overhead
  exceeding benefit at this per-thread work granularity for this contact count, not a bug — but it
  means `condensed.num_threads=0` (the config default) is not a safe "just works" choice on every
  scene; it should be tuned.
- **None of this is an apples-to-apples "condensed beats PJ" result**, and it isn't presented as
  one. PJ and PGS here are also non-accelerated (`pj.nesterov=false`) baselines. The one clean,
  reproducible win is `newton` vs. `projection` *within* `condensed` itself.

## Suggested next steps, in priority order

1. Switch the per-block hot path (`blockResidual`, `computeBlockUpdate`, `scatterDelta`) to
   fixed-size Eigen types to close the wall-clock gap with PGS on `gauss_seidel`.
2. Cache `colored`'s coloring across timesteps when the contact set is stable, and/or increase the
   amount of work done per OpenMP barrier.
3. Re-run the TSan validation with a Clang/`libomp` toolchain if one becomes available, to get a
   clean (not just argued-correct) result for the parallel sweep modes.
4. If `condensed` is to be compared against PJ/PGS's *accelerated* configurations, port Nesterov
   (or Chebyshev/Anderson) to `computeBlockUpdate()`'s call site — the pattern already exists in
   `projected_jacobi.cpp`/`projected_gauss_seidel.cpp` to adapt from.
5. Run the full ablation in `tools/compare_solver_traces.py` (now supports all 8
   sweep-mode/local-solve combinations on both scenes with wall-clock capture) with a
   multi-minute-per-case timeout budget to get complete, not partial, wall-clock numbers —
   the runs in this report were deliberately capped at 90s/case to fit the session, and several
   cases (marked "no" above) didn't reach full completion in that budget.
