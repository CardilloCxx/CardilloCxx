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

**This was true when first asked, and is the starting point for the rest of this report** (the
"Benchmark results" section below predates the change described next): all four sweep modes were
plain fixed-point iterations, with relaxation (`pj_relaxation`) and per-contact step damping
(`pj_alpha`) as the only step-size controls, identical to what PGS already had, and no Nesterov/
Chebyshev/Anderson.

**As of the "Update" section below, Nesterov acceleration is implemented** for `gauss_seidel`,
`jacobi`, and `colored` (not `chaotic`, deliberately — see that section), reusing PJ/PGS's existing
`pj.nesterov`/`pj.nesterov_beta_threshold`/`pj.nesterov_restart_limit` config keys rather than
adding new ones. Chebyshev and Anderson are still not ported; see "Suggested next steps."

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

## Update: Nesterov acceleration + fixed-size buffer optimization

Follow-up work, done in this order because the user specifically asked to check whether Nesterov
(PJ's biggest documented win) could be ported before anything else.

### Nesterov acceleration, generalized across sweep modes

Added a generic Nesterov (FISTA-style) outer loop in `CondensedSolver::solve()`, ported
line-for-line from `ProjectedGaussSeidelSolver`'s own Nesterov branch (momentum extrapolation
between sweeps, adaptive restart on residual growth / excessive momentum coefficient / direction
reversal, momentum permanently disabled after `pj.nesterov_restart_limit` restarts). It's
implemented once, generically, via a `doSweep(lambda, u_corr)` closure that dispatches to whichever
of `gaussSeidelSweep`/`jacobiSweep`/`coloredSweep` is configured — Nesterov only ever touches the
outer state between sweep calls, never the sweep internals, so this required no changes to the
sweep functions themselves. No new config keys: it reuses `pj.nesterov`,
`pj.nesterov_beta_threshold`, `pj.nesterov_restart_limit`, the same keys PJ and PGS already read.

`chaotic` is deliberately excluded — its state lives in an atomic array, not a plain vector, and
momentum-extrapolating a value that's already being intentionally kept stale by design is a
separate question not tackled here.

**First bug caught, fixed immediately:** the initial implementation reassigned `lambda`/`u_corr`
from the lagging `lambda_k`/`u_k` variables after the loop exited. PGS's own code doesn't do this —
`lambda`/`u_corr` already hold the correct latest post-sweep values from inside the loop (the sweep
functions mutate them by reference), and `lambda_k`/`u_k` are one iteration behind at that point.
Caught by re-reading the ported reference code line-by-line before running it, not by test failure.

**Correctness check** (domino, `gauss_seidel`, `pj.nesterov=true` vs `false`): the accelerated
result (`totalKE=0.814479477909828`) matches PGS's own Nesterov-accelerated result
(`0.814479477909845`, both scenes inherit `pj.nesterov=true` from `scene_standard.config`'s
default) to 11 significant digits, and the non-accelerated result reproduces the pre-Nesterov
number exactly (`0.811464242457628`, bit-identical to before this change — a clean regression
check). Sweep count: 11979 (plain) → 3057 (Nesterov) on domino `gauss_seidel`, landing on
essentially the same iteration count as PGS's own Nesterov run (3057) — expected, since the
underlying per-sweep math is identical between the two. `jacobi`+Nesterov similarly converges to
PJ's own fixed point (`totalKE=0.814196987085242` vs PJ's `0.814196987077336`) in 3057 sweeps vs
PJ's 3056.

**A real, scene-dependent limitation found — and then root-caused, not just worked around:**
`jacobi`+Nesterov originally diverged on slinky (`[Condensed] Divergence detected after 949
iterations. Residual norm: -nan`), while `gauss_seidel`+Nesterov was stable there. The first
write-up of this blamed Nesterov's momentum shrinking Jacobi's already-smaller stability margin —
that explanation was wrong. Testing a **plain, unaccelerated** `jacobi` sweep (Nesterov off
entirely) on slinky at the same `pj.alpha=0.3` showed it diverges too, identically. Momentum was
never the cause: `slinky`'s default `pj.alpha=0.3` is tuned for `gauss_seidel`/`colored`, and
`jacobi` has a smaller stability margin at a given alpha than those do — the exact same
per-sweep-mode instability already documented in `domino`'s own `scene_condensed.config` (which
needed `pj.alpha=0.001`, not the header default of 0.3, for `jacobi` to be stable there). Dropping
`pj.alpha` to `0.02` made `jacobi`+Nesterov stable on slinky over 1000+ steps (sustained ~11 it/s,
zero divergence) — the fix is retuning alpha per sweep mode, not disabling Nesterov.

Separately, while investigating this, a genuine bug was found and fixed in the Nesterov loop
itself: the residual check driving the loop's exit/restart decision was routed through a helper
that throws `std::runtime_error` immediately on a `nan`/`inf` residual. That made the loop's own
`if (!std::isfinite(err)) restart = true;` branch permanently dead code — the exception fired
before the restart logic could ever run, even though that logic (ported from PJ's own
`nesterov_loop()`, which never throws) exists specifically to recover from exactly this case. Fixed
by computing the residual without throwing inside the Nesterov loop, letting restart/momentum-
disable react to a non-finite value exactly as it reacts to a merely-growing one, and only raising
(matching PJ/PGS's existing, inherited behavior) once momentum is fully disabled and the plain
sweep output is itself still non-finite — i.e. once there is genuinely nothing left to restart to.
Also added a finiteness check on the freshly-extrapolated `(lambda_y, u_corr_y)` state before it is
fed into the next sweep: `betak1` is bounded to `[0,1]`, but the extrapolated vector's magnitude is
not, so an overshoot can still produce a non-finite state even with a well-behaved coefficient.
**Practical guidance updated: retune `pj.alpha` down when switching a scene's sweep mode to
`jacobi` (accelerated or not); the divergence is a per-sweep-mode stability property of the
scene/alpha combination, not a Nesterov-specific risk.**

### Fixed-size buffer optimization (the top item from the original "next steps" list)

Replaced every per-block temporary in the hot path (`blockResidual`'s residual, `lamOld`, the
Newton/projection update, the scatter delta, the Jacobi-gather accumulator, the chaotic-mode
per-block scratch) with `Vectorr<6>` — a fixed-size, stack-allocated Eigen type already used
elsewhere in this codebase — instead of dynamically-sized `VectorXr`. `RowBlock::Ja`/`Jb`/`Gii`/
`GiiInv` were deliberately left as dynamically-sized `MatrixXXr`: they're built once per `solve()`
call in `CondensedAssembler`, not once per block per sweep, so they were never the bottleneck.
Every block's row count is ≤ 6 (1 for frictionless contacts, 3 for frictional, up to 6 for
springs/dampers) and every body's DOF count is ≤ 6, so a fixed 6-element buffer with `.head(n)`
views covers every case exactly. `projectBlock()` was templated on `Eigen::MatrixBase<Derived>` so
the same function body works on both `VectorXr::segment()` views and `Buf6::head()` views.

**Correctness check**: re-ran `gauss_seidel`/`jacobi`/`colored` on domino at the exact same
tolerance as before the refactor — all three reproduced their pre-refactor kinetic energy and
sweep count **bit-for-bit** (e.g. `gauss_seidel`: `0.811464242457628`, 11979 sweeps, identical to
15 displayed digits). `chaotic` landed on a nearby but not identical value (expected — it's
randomized), with a comparable sweep count.

**Performance result**, domino, `pj.nesterov=true`, sequential (non-contaminated) runs, 3 steps:

| Case | Wall-clock before | Wall-clock after | Sweeps |
|---|---|---|---|
| PJ (unmodified, for reference) | 6.76s | 6.76s (unchanged) | 3791 |
| Condensed `gauss_seidel`+Nesterov | 19.16s | 13.64s (−29%) | 4149 |
| Condensed `jacobi`+Nesterov, 4 threads | 8.37s | 5.18s (−38%) | 3796 |

A thread-count sweep for `jacobi`+Nesterov post-optimization (sequential runs, no cross-run
contention) found 8 threads as the sweet spot on this 16-core machine:

| Threads | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| Wall-clock | 11.11s | 5.95s | 3.55s | **3.02s** | 4.08s |

**PJ is no longer the fastest option on domino**: PJ takes 5.35s on a clean sequential run; condensed
`jacobi`+Nesterov at 8 threads takes 3.02s — **~1.8x faster than PJ**. On slinky, `gauss_seidel`+
Nesterov (sequential, no threading) takes 8.39s vs. PJ's 12.04s — **~1.4x faster than PJ**, with
`jacobi`+Nesterov unusable there (see the divergence finding above). Both wins are with
`condensed.local_solve=projection`; `newton` on top of an already-Nesterov-accelerated sweep did
not help further in this round of testing (slightly slower than projection alone on slinky's
already-fast-converging accelerated run — the per-contact Newton cost isn't repaid when only ~30
sweeps are needed in total).

**Important caveat on these numbers**: they're from short runs (3 steps on domino, 200 steps on
slinky) picked to fit within the session, on one machine, with one thread count tried exhaustively
(domino) and one (8) assumed-reasonable for slinky. They demonstrate that condensed *can* beat PJ,
not that it always will — a proper claim would need the full `tools/compare_solver_traces.py`
ablation (now updated with `condensed` cases) run to completion across both scenes with a
generous timeout, which is listed below as still-open work.

## Update: large `pj.alpha` for `gauss_seidel`/`colored` (not `jacobi`)

Separate from Nesterov, the user found empirically that `condensed.sweep_mode=gauss_seidel` and
`colored` tolerate — and benefit from — much larger `pj.alpha` than `jacobi` does at the same scene
and `moreau.theta`. On `domino` (`colored`, 8 threads), `pj.alpha` swept from 0.1 to 1.0 monotonically
*improved* wall-clock (0.1: 131s, 0.25: 119s, 0.5: 117s, 1.0: 101s), the opposite of `jacobi`'s own
sensitivity on the same scene (0.05: 91s, 0.1: 97s — smaller is better there). Confirmed on today's
code: `condensed` `colored`+`alpha=1.0`+8 threads takes ~98.4s on `domino`'s full `scene_condensed.config`
run vs PJ's ~217.6s at PJ's own tuned `alpha=0.1` (**~2.2x**) — and on `slinky`, PJ at its own tuned
`alpha=0.3` didn't even converge within a 900s budget in one config tested (its Nesterov loop's debug
trace showed a residual barely moving, then trending upward, over tens of thousands of iterations),
while `condensed` `colored`+`alpha=1.0` completed the full 1000-step run in ~106.8s. **PJ itself is
*not* tolerant of large alpha** — a PJ run at `alpha=1.0` on `domino` didn't finish 300 steps within a
10-minute budget, confirming this is specific to `condensed`'s `gauss_seidel`/`colored` sweep
structure, not a general property of the underlying iteration.

Caveat: pushing `jacobi` to match `gauss_seidel`/`colored`'s large-alpha tolerance does not work —
`jacobi` needs a *smaller* alpha than the scene's `gauss_seidel`/`colored`-tuned default (see the
Nesterov update above for the concrete `slinky` numbers), consistent with Jacobi's smaller stability
margin at a given alpha. Domino's and slinky's `scene_condensed.config` reference configs were
updated to `pj.alpha=1.0` with `condensed.sweep_mode=colored` as their shipped defaults.

## Update: true Schur-complement elimination of the compliant chain (`condensed.true_schur`)

### Motivation

Investigating why PJ sometimes needed dramatically fewer sweeps than condensed on `slinky` (and
separately, why a specific `slinky` config drove PJ's own Nesterov loop to churn for 60,000+
iterations without converging — see the earlier "improve this further" discussion) pointed at a
structural gap: every condensed sweep mode treats every row — spring, damper, frictionless
contact, frictional contact — as an independent block, and lets the *outer sweep loop itself*
diffuse information between rows that share a body. For `slinky`'s 825-segment chain (~824
spring/damper rows forming a path graph), that diffusion — not the frictional contact
nonlinearity — was the dominant cost. PJ doesn't have this problem because it factorizes and
re-solves the *entire* system every sweep, seeing the chain's exact coupling every time, at the
cost of a much more expensive per-sweep operation (a full triangular solve, not a condensed
update).

The user had tried one thing here already and confirmed (after I searched the repo and couldn't
find it) that it was narrowly: swapping `condensed_assembler.cpp`'s `GiiInv` from the diagonal-only
inverse to the full local block inverse (`blk.Gii + diag(complianceDiag)`, inverted via
`partialPivLu`) — code that still exists as a commented-out TODO in `updateCompliance()`. That is
**not** a Schur complement — it only removes coupling *within* one row's own components (e.g. a
frictional contact's normal/tangential coupling), not the coupling *between* different rows that
share a body, which is the actual chain-diffusion bottleneck. It's also the exact bug already
documented above (2.8x kinetic-energy error on domino from using the full inverse instead of the
diagonal one matching `PgsAssembler::DinvDiag()`) — explaining why that attempt "didn't seem
correct."

### The actual construction

Working from the same `S = W·Minv·W^T + diag(C)` system every PGS-family solver already solves,
partition rows into bilateral (springs+dampers, pure equalities — `projectBlock()` already no-ops
for these) and contact (subject to the friction-cone projection). The Schur complement over
contacts alone is `(Scc − Scb·Sbb⁻¹·Sbc)·λc = rhsc − Scb·Sbb⁻¹·rhsb`. Explicitly forming this
(generally dense) matrix is exactly what this solver family exists to avoid, and isn't necessary:
because bilateral rows are pure equalities, whatever the current contact impulses are, there is a
single exact `λb` that zeroes every bilateral row's residual at once — `δλb = Sbb⁻¹·r_old`, where
`r_old` is the *already-existing* `blockResidual()` evaluated at the bilateral blocks' current
state. Applied via the *already-existing* `scatterDelta()`, this reuses every existing primitive —
no new sign convention, no rhs()/blockResidual()/scatterDelta() changes. Each outer iteration
becomes: (1) exact bilateral step (one `Sbb` solve via a precomputed factorization), then (2) the
usual sweep restricted to just the contact blocks.

New, generic, reusable utility: `BlockSparseLDLT` (`src/misc/block_sparse_ldlt.{hpp,cpp}`) — a
block-sparse LDLT factorization over a graph of small dense blocks, never a dense/global matrix.
Verified in isolation (a standalone test compiled outside the main build, comparing against a known
solution and against Eigen's own dense LDLT) on a chain, a branching hub-and-leaves graph in two
different elimination orders (exercising real fill-in), and a 270-node stress graph mirroring
`hangbridge`'s rope+plank topology — all four to machine precision (~1e-16 relative error).

**A real bug found and fixed during this work, not just during isolated unit testing**: the first
integration attempt used natural (row creation) order for elimination — correct in principle, and
zero-fill-in on every pure-chain scene — but it **hung** on `hangbridge`, the one example scene
with a real branching compliant network (tripod apexes / deck planks, degree 4-6; confirmed via an
Explore-agent survey of every example scene's spring/damper topology). Natural order is not
fill-reducing for a branching graph, and the fill-in cascade made factorization impractically slow.
Fixed by computing elimination order internally via greedy minimum-degree (a standard, well-known
sparse-Cholesky ordering heuristic) instead of trusting the caller's order — it recovers exactly
the chain order (zero fill-in) on every path-graph scene, so nothing is given up on the common case,
and the 270-node hangbridge-shaped stress test now factors in under a millisecond.

### Validation

- **Domino** (no springs): `condensed.true_schur=true` vs `false` — bit-identical `totalKE`
  (`8.96934551302463`) and position fingerprint. Required by construction (nothing to eliminate
  when there are zero bilateral rows) and confirmed empirically, the critical regression guardrail.
  This holds up unchanged.
- **Slinky — corrected after the user asked for a re-test that found the same order-of-magnitude
  with the feature on or off.** The first pass through this section claimed a "~25.8x speedup,
  ~2000/step to ~2/step" finding. That was **wrong** — an artifact of reading a truncated log tail
  that happened to show only the first simulated step. Re-measured properly (full per-step logs,
  three independent re-runs, and a from-scratch build of the pre-`true_schur` commit to rule out an
  unrelated regression as the explanation): `slinky`'s coiled geometry means the very first step has
  **zero** contacts (pure spring-chain relaxation) — there, `true_schur` genuinely needs 2 outer
  iterations vs 2532 without it, confirming the theory exactly where it applies. But contacts start
  forming by the *second* step (414) and grow rapidly (504 → 786 → 1005 → 1206 → 1506 → 2154 → 2532
  within the first 9 steps, as the coil self-collapses). From that point on, `true_schur` on vs off
  gives essentially the same iteration count per step (e.g. step 7: 16393 vs 16243; step 13: 46417
  vs 46328) — contact resolution, which `true_schur` does not touch, dominates almost immediately.
  Net effect over a full 20-step run: **109.8s with `true_schur=true` vs 85.4s without — a
  slowdown**, because the per-iteration `LDLT` solve cost is paid every step but only repaid on the
  one step that's actually contact-free.
- **Hangbridge** (branching topology, contact-heavy from the very first step: 232 contacts on 441
  bodies): correctness verified — `totalKE` agrees with a direct PGS run to ~1e-6 relative both with
  `true_schur=false` (bit-identical to PGS: `0.00350788969155182`) and `true_schur=true`
  (`0.00350788582097124`). Performance was roughly a wash to slightly worse (0.97s vs 0.57s for 20
  steps) — consistent with the corrected slinky finding above, not a separate anomaly: this scene
  never has a contact-free step to repay `true_schur`'s per-iteration cost.

**Practical guidance, corrected**: `condensed.true_schur=true` is a real win specifically for
**contact-free** portions of a simulation (a compliant structure settling before it touches
anything, or genuinely between contact events) — not a general win for "scenes with a compliant
chain," which was the original (too broad) framing. For a scene that's in contact for most of its
runtime, including `slinky` almost immediately, expect it to be neutral to a net slowdown, not a
dramatic speedup. Default is `false`. Before trusting a number on a new scene: look at iteration
counts *per step*, not a run's total or its first logged value.

**How the original mismeasurement happened, for the record**: the original benchmark command
piped output through `tail -6`/similar to keep the transcript short, which — for a run whose early
steps print short lines and later steps print long wrapped ones — silently showed only the first
step's numbers, not a representative sample. The fix going forward: when characterizing
per-step behavior, always inspect the full step-by-step sequence (e.g. `grep -o "iters=[0-9]*"`
over the whole log), never a tail/head slice, especially for a scene whose difficulty is expected
to change over the run (contacts forming, a structure settling, etc.).

## Update: fixed-size BlockSparseLDLT internals + cached elimination order

Follow-up work after `wilberforce` (a pure DAE, `collision.disable_all=true`) showed `true_schur`
"brilliant" there — user request to push further: compile-time-bounded block dimensions, and
symbolic/numeric caching.

**Fixed-size blocks**: tried in two places. `BlockSparseLDLT`'s own internal storage (used by
`factor()`/`solve()`, called once per pivot/once per call) switched from heap-allocated `MatrixXXr`
to fixed `Matrixr<6,6>` buffers accessed via `.topLeftCorner()` — verified correct (standalone
stress test, ~1e-16 accuracy) and a genuine ~4% wall-clock win on `wilberforce` (bilateral blocks
are dim=6, a good fit). Applying the *same* treatment to `RowBlock` itself (`Ja`/`Jb`/`Gii`/
`GiiInv`, read by every sweep function) was also tried — measured a **~20% wall-clock regression on
`domino`** (own from-scratch build vs an A/B-tested stash of the pre-change commit, not a rough
estimate). Root cause: domino has ~13k mostly small-dim (1 or 3) contact blocks; forcing them into
a 6x6 buffer wastes most of that memory as unused padding, and reading that padding-heavy struct
thousands of times per sweep hurts cache locality more than avoiding heap allocation helps. Reverted
that part; kept the `BlockSparseLDLT`-only version, which doesn't have this problem (bilateral
blocks tend to actually use most of a 6x6 buffer).

**Cached elimination order**: `CondensedAssembler` now caches the minimum-degree order across
`solve()` calls, keyed on the bilateral graph's structure (dims + coupled-block pairs), invalidated
only if that structure changes (never does, in every current scene — constraints aren't created or
destroyed at runtime). `factorWithOrder()` still does the full numeric factorization on the current
step's values regardless of a cache hit, so this can only affect performance, never correctness —
verified bit-identical `totalKE` on domino/slinky/wilberforce/hangbridge before and after. Measured
~14% faster `Condensed Setup` on `wilberforce` (327.9µs → 282.6µs average over 2000 calls). Total
wall-clock effect was much smaller, because `Condensed Setup` itself is only ~5% of that scene's
runtime — `Output Write` (VTK I/O, `output.interval_steps=1` writing every step) is ~90%. Noting
this plainly rather than folding it into a bigger-sounding total-wall-clock number: if wall-clock on
a similar scene matters, output frequency is a bigger lever than anything in this solver.

**A real correctness question raised during review, resolved with a concrete test, not just an
argument**: does this cache actually invalidate correctly if the bilateral topology changes at
*runtime* (a constraint added/removed mid-simulation via `World::markStructureDirty()`), given no
existing example scene's topology ever changes? Traced the mechanism:
`DynamicsAssembler::updateStateDependentTerms()` calls `rebuildInteractionW_()` (which fully rebuilds
`constraintResults()`) unconditionally every step, regardless of `m_structure_dirty` — that flag only
gates the separate, cheaper-to-skip body-offset/dof bookkeeping in `refreshState()`. So
`CondensedAssembler::buildTopology()` always reflects the live constraint set, and the cache's own
`dims == cached && edgeNodes == cached` check (computed fresh from that live data every call) is a
sufficient, self-contained invalidation signal — it doesn't need to consult `m_structure_dirty` at
all. Rather than resting on that argument, added `examples/scenes/dynamic_constraint`: three point
masses, one spring from the start, a second added via `updateScene()` at t=0.02s. Confirmed via
`debug.pj=true` that `nSprings` actually grows 1→2 mid-run, and that `condensed` with
`true_schur=true` (cache exercised) produces a `totalKE`/`posNormSum` fingerprint identical to both
PGS and `condensed` with `true_schur=false` across the structural change. A genuine regression test
for this guarantee, not just a plausibility argument.

## Suggested next steps, in priority order

1. ~~Harden Nesterov against the `jacobi` divergence case~~ — done: root-caused (per-sweep-mode
   alpha stability, not a momentum effect) and the dead-restart-branch bug that made the loop
   unable to even attempt recovery is fixed. See the update above.
2. Run the full ablation in `tools/compare_solver_traces.py` with a multi-minute-per-case timeout
   to replace this report's short, hand-picked comparisons with complete, systematic ones —
   including a proper per-scene thread-count sweep (only domino's was swept exhaustively above),
   and a `jacobi`+Nesterov slinky run at the corrected `pj.alpha=0.02`.
3. ~~Cache `colored`'s coloring across timesteps~~ — turned out not to matter: at the large-alpha
   operating point, `CondensedColoring` is only ~1.2% of wall-clock (iteration counts collapsed
   from tens of thousands to double/triple digits once alpha was tuned properly), so the bottleneck
   moved to the sweep itself before this was ever implemented. Not worth doing now.
4. Re-run the TSan validation with a Clang/`libomp` toolchain if one becomes available, to get a
   clean (not just argued-correct) result for the parallel sweep modes.
5. Consider Chebyshev/Anderson acceleration too (PJ has both) now that the generic `doSweep()`
   wrapper pattern exists for bolting an outer accelerator onto any sweep mode.
6. Investigate why `true_schur` doesn't stack with the large-alpha/colored/Nesterov tuning that
   already works well — is there a way to get both the exact-elimination win *and* the
   large-alpha/parallel win on the same run, or are they fundamentally solving the same problem
   twice? Would need to understand what's actually still costing hundreds of iterations at
   `alpha=1.0` if it isn't chain diffusion.
7. Extend `true_schur`'s exact elimination to cover contact-contact coupling *through* the chain
   too (a fuller static-condensation/Guyan-reduction treatment), not just bilateral-bilateral
   coupling — may be what's needed to make it help on contact-heavy branching scenes like
   `hangbridge` where iteration count is dominated by something `true_schur` today doesn't touch.
8. **Compile-time block dimensions** (flagged by the user as a distinct future step, not started):
   extend the `Buf6`/`Vectorr<6>` fixed-size-buffer optimization already applied to per-block
   *vectors* in `condensed_solver.cpp`'s hot path to `RowBlock`'s *matrices* (`Ja`, `Jb`, `Gii`,
   `GiiInv`), which are still dynamically-sized `MatrixXXr` built once per `solve()` call. Lower
   risk to attempt after `true_schur` is further validated, not alongside it — bundling two
   structural changes at once multiplies the surface area for a subtle bug, which is exactly what
   burned the user's own single-line `GiiInv` attempt earlier in this document.
