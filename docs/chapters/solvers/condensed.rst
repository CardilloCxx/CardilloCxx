Condensed (Block-Sparse) Solver
================================

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

``condensed`` is a family of contact solvers sharing one assembler and one
solver class, selected entirely through config keys rather than separate
backends. It generalizes :doc:`projected_gauss_seidel` (PGS) along two
independent axes:

- **sweep strategy** -- how the outer iteration visits contact/constraint
  blocks: sequential (``gauss_seidel``), embarrassingly parallel
  (``jacobi``), graph-colored parallel (``colored``), or deliberately
  unsynchronized parallel (``chaotic``);
- **local solve** -- how one frictional contact block is updated within a
  sweep: the usual linear-step-then-project (``projection``), or a guarded
  semismooth Newton solve (``newton``) that can converge a single contact to
  machine precision in a handful of iterations instead of tens.

Motivation: unlike :doc:`projected_jacobi` (PJ), which factorizes the full
mass+bilateral system once but re-solves it (a global sparse triangular
back-substitution) on **every** sweep, ``condensed`` -- like PGS -- never
touches a global system matrix after setup. Unlike PGS, it is not hardwired
to one sequential order or one local-solve rule, so the sweep/local-solve
choice can be tuned per scene, and the sweep loop can use OpenMP.

**Acceleration:** ``gauss_seidel``, ``jacobi``, and ``colored`` support Nesterov
acceleration via the same ``pj.nesterov`` config key PJ and PGS already use
(no new keys) -- see `Nesterov acceleration`_ below, including a real
scene-dependent divergence case to be aware of. Chebyshev and Anderson (both
PJ-only today) are not yet ported.

Block-sparse representation
----------------------------

The assembler (``CondensedAssembler`` / ``RowBlock`` / ``CondensedTopology``
in ``src/physics/assembly/condensed_assembler.hpp``) never builds an
``Eigen::SparseMatrix``, not even for matrix-free matvecs. Every spring,
damper, frictionless contact (1 row), or frictional contact (3 rows: normal +
2 tangential) becomes a ``RowBlock``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Meaning
   * - ``Ja``, ``Jb``
     - Dense Jacobian rows for this block's two sides (``dim`` x ``aDof``/``bDof``), zero-sized if that side is static
   * - ``Gii``
     - Local (mass-only) Delassus block, :math:`G_i = J_a\,\mathrm{diag}(M^{-1}_a)\,J_a^\top + J_b\,\mathrm{diag}(M^{-1}_b)\,J_b^\top`
   * - ``complianceDiag``
     - :math:`C_\mathrm{row}/(\theta\,dt^2)` (springs), :math:`A_\mathrm{row}/(\theta\,dt)` (dampers), zero (contacts)
   * - ``GiiInv``
     - Diagonal-only inverse of :math:`G_i + \mathrm{diag}(\mathrm{complianceDiag})` -- see note below
   * - ``bodyIndexA/B``, ``aOff/bOff``, ``aDof/bDof``
     - Which two bodies this block touches and their velocity-vector offsets (``-1``/``0`` if static)

.. note::
   ``GiiInv`` is a **diagonal**-only approximation of :math:`G_i^{-1}`, matching
   ``PgsAssembler::DinvDiag()`` -- the function PGS's own sweep loop actually
   calls at runtime. ``PgsAssembler::Dinv()`` (the full dense block inverse)
   is dead code for contacts in this codebase; only the bilateral-only
   :doc:`conjugate_gradient` solver uses it. Using the full block inverse
   here instead of the diagonal one changes the converged answer, not just
   the convergence rate -- confirmed by direct comparison against PGS during
   development (a divergent-looking-but-"converged" 2.8x kinetic energy
   error on the domino scene). ``Gii`` itself is still stored as the full
   3x3/6x6 matrix because the semismooth Newton local solve needs the true
   local coupling; only the *linear* local-solve path uses the diagonal
   approximation.

Two body-index maps fall out of the same block list for free: a per-body
list of incident blocks (``CondensedTopology::blocksOfBody``, used by
``jacobi``'s gather pass and to build ``colored``'s adjacency graph) and the
packed row layout (springs, then dampers, then frictionless contacts, then
frictional contacts -- identical to PGS's own layout, reusing
``DynamicsAssembler``'s own contact row ordering (``impulse_base_index``)
rather than recomputing it).

Sweep modes
-----------

All four modes share the same per-block update
(:math:`\lambda^{(k+1)} = \mathrm{proj}\bigl(\lambda^{(k)} + \omega\,\alpha\,D_i^{-1}\,r_i\bigr)`,
:math:`r_i` the block's local residual) -- they differ only in *when* each
block reads/writes the shared velocity-correction state ``u_corr``.

.. mermaid::

   flowchart TD
       A[Build RowBlock topology + compliance] --> B[Warmstart lambda, gather initial u_corr]
       B --> C{sweep_mode}
       C -- gauss_seidel --> D[Sequential: propagate u_corr immediately after each block]
       C -- jacobi --> E[Pass 1: update all blocks vs frozen u_corr snapshot -- Pass 2: each body gathers from its own incident blocks]
       C -- colored --> F[Greedy-color the shared-body adjacency graph once -- parallel within a color, sequential across colors]
       C -- chaotic --> G[Shuffle block order periodically -- parallel over the permutation, atomic u_corr accumulator]
       D --> H{Converged?}
       E --> H
       F --> H
       G --> H
       H -- No --> C
       H -- Yes --> I[Store impulses for warmstart, return u_free - u_corr]

**gauss_seidel** (default)
   Sequential, in the blocks' natural (packed) order, propagating ``u_corr``
   immediately -- structurally identical to :doc:`projected_gauss_seidel`'s
   own loop, just expressed over ``RowBlock`` s instead of a
   ``BlockDiagonal``/sparse ``W`` pair. Never parallelized (there is nothing
   embarrassingly parallel about a sequential update by construction).

**jacobi**
   Two OpenMP passes per sweep, run inside a single ``#pragma omp parallel``
   region (one thread-team spawn per sweep, not two): pass 1 updates every
   block's own impulse against a *frozen* snapshot of ``u_corr`` from the
   previous sweep, writing only to that block's own (disjoint) slice of
   ``lambda``/``dlambda``; pass 2 has each *body* gather corrections from
   only its own incident blocks (via ``blocksOfBody``) into its own slice of
   the new ``u_corr`` -- a gather, not a scatter, so this needs no atomics at
   all, not even for the reduction.

**colored**
   Greedy Welsh-Powell coloring (``src/misc/graph_coloring.hpp``) of the
   *shared-body* adjacency graph (two blocks are adjacent iff they touch the
   same body) -- deliberately **not** the dense Schur-complement/Delassus
   sparsity pattern, which for a compliant chain would be fully dense and
   collapse every color into one. Blocks within a color are updated in
   parallel (no two share a body, so no two write the same ``u_corr``
   entries); colors are processed sequentially, one ``#pragma omp for``
   each, inside a single ``#pragma omp parallel`` region for the whole
   sweep.

**chaotic**
   The experimental mode. Block processing order is a permutation, reshuffled
   every ``condensed.chaotic_reshuffle_interval`` sweeps; the permutation
   property means each block index is touched by exactly one thread per
   sweep, so writes to ``lambda`` stay race-free. ``u_corr`` is different: two
   blocks sharing a body genuinely can update it concurrently, and that
   staleness is the intended effect (Chazan-Miranker chaotic relaxation, not
   a bug), implemented with a real
   ``std::vector<std::atomic<real_t>>`` (relaxed loads/``fetch_add``) --
   never a plain-``double`` data race, which would be undefined behavior in
   C++, not just imprecise.

.. warning::
   ThreadSanitizer flags the OpenMP regions in ``jacobi``/``colored``/
   ``chaotic``. This was tracked down to a documented, decade-old limitation:
   GCC's ``libgomp`` is not itself TSan-instrumented, so TSan cannot see its
   barrier/team-join synchronization and reports the (correctly
   barrier-ordered) accesses on either side as a race. A minimal standalone
   reproducer (two plain ``#pragma omp for`` loops in one
   ``#pragma omp parallel``, disjoint writes, no other synchronization
   primitive involved) triggers the identical report under this
   GCC/``libgomp``/TSan combination, confirming this is not specific to this
   solver's code. No Clang/``libomp`` toolchain was available in the
   development environment to get a cleanly-instrumented run. Correctness
   confidence instead rests on the disjoint-access arguments documented
   inline in ``condensed_solver.cpp`` plus bit-identical results across 1,
   4, and 16 threads on both benchmark scenes.

Local solve: projection vs. semismooth Newton
-----------------------------------------------

**projection** (default) clips exactly as :doc:`projected_gauss_seidel`
does: normal impulse clamped to :math:`\Lambda_N \le 0`, tangential impulses
scaled toward the Coulomb disk of radius :math:`-\mu\,\Lambda_N`.

**newton** replaces that clip, for frictional (3-row) blocks only, with a
guarded Alart-Curnier semismooth Newton solve
(``src/physics/solver/local_contact_newton.{hpp,cpp}``). Given the block's
local Delassus matrix :math:`G_i` and current residual :math:`r`:

.. math::

   \rho_N = 1/G_{i,00},
   \qquad
   \rho_T = 1/\lambda_\mathrm{max}(G_{i,TT})

(:math:`G_{i,TT}` the 2x2 tangential sub-block; :math:`\lambda_\mathrm{max}`
via the closed-form 2x2 eigenvalue, no ``SelfAdjointEigenSolver`` needed).
With :math:`M = I - \mathrm{diag}(\rho_N,\rho_T,\rho_T)\,G_i` and
:math:`c = \mathrm{diag}(\rho)\,(r + G_i\,\lambda)`, the trial state
:math:`y(\lambda') = M\lambda' + c` is case-split into separation
(:math:`y_N \ge 0 \Rightarrow g = 0`), sticking
(:math:`\|y_T\| \le -\mu y_N \Rightarrow g = y`), or sliding
(:math:`g_T = -\mu y_N\,\hat t`, :math:`\hat t = y_T/\|y_T\|`) to build
:math:`\Phi(\lambda) = \lambda - g(\lambda)` and its Jacobian, then takes a
Newton step :math:`J_\Phi\,\delta = -\Phi`.

The solve is **guarded**: it reports failure (never leaves ``lambda`` in a
partially-updated state) on a singular/near-singular local block, a
non-finite step, a step that increases :math:`\|\Phi\|` after the first
iteration, or exhausting ``condensed.newton_max_iters`` without meeting
``condensed.newton_tol``. On failure the caller falls straight through to
the ``projection`` update for that block -- ``newton`` never risks being
*less* robust than ``projection``, only sometimes not worth its extra cost.

True Schur-complement elimination of the compliant chain (``condensed.true_schur``)
-------------------------------------------------------------------------------------

**The problem this solves.** Every sweep mode above treats *every* row --
spring, damper, frictionless contact, frictional contact -- as an
independent ``RowBlock``, and lets the outer sweep loop implicitly diffuse
information between rows that share a body. For a scene with a long
compliant chain (``slinky``: 825 segments, ~824 spring/damper rows forming a
path graph), that diffusion is the dominant cost: a contact's effect on one
end of the chain has to slowly propagate, one Gauss-Seidel row-visit at a
time, before a contact at the other end can "feel" it. :doc:`projected_jacobi`
does not have this problem -- it factorizes and re-solves the *entire*
system every sweep, so it sees the chain's exact coupling on every
iteration -- but pays for that with a much more expensive per-sweep cost
(a full triangular back-substitution, not a condensed block update).
``condensed.true_schur=true`` gets ``condensed`` PJ's exactness on the
bilateral (spring+damper) rows without paying PJ's per-sweep cost, by
eliminating those rows in closed form once per outer iteration instead of
letting them participate in the sweep at all.

**The math**, worked from the same ``S = W*Minv*W^T + diag(C)`` system every
PGS-family solver already solves (see ``PgsAssembler::rhs()``/``Dinv()`` and
their condensed transcriptions -- nothing below is re-derived independently
of that verified ground truth). Partition rows into bilateral (springs +
dampers, pure equality constraints -- ``projectBlock()`` already no-ops for
these) and contact (subject to the friction-cone projection):

.. math::

   S = \begin{bmatrix} S_{bb} & S_{bc} \\ S_{cb} & S_{cc} \end{bmatrix},
   \quad
   \lambda = \begin{bmatrix} \lambda_b \\ \lambda_c \end{bmatrix}

Eliminating :math:`\lambda_b` gives the Schur complement over contacts alone:
:math:`(S_{cc} - S_{cb} S_{bb}^{-1} S_{bc})\,\lambda_c = \text{rhs}_c - S_{cb} S_{bb}^{-1}\,\text{rhs}_b`.
Explicitly forming this (generally dense, since :math:`S_{bb}^{-1}` is dense
even though :math:`S_{bb}` is sparse) matrix is exactly what this solver
family exists to avoid, and turns out not to be necessary: because bilateral
rows are pure equalities with no projection, whatever the current contact
impulses are, there is a single exact :math:`\lambda_b` that zeroes every
bilateral row's residual simultaneously --
:math:`\delta\lambda_b = S_{bb}^{-1}\,r_{\text{old}}`, where :math:`r_{\text{old}}`
is the *already-existing* ``blockResidual()`` evaluated at the bilateral
blocks' current state. Applied via the *already-existing* ``scatterDelta()``,
this is architecturally just "one generalized block update over the whole
bilateral row-set as a single joint block, using its exact factorized
inverse" -- the same delta-then-scatter pattern every sweep function already
uses per-row, just applied jointly instead of per-row-and-iteratively. No
rhs()/blockResidual()/scatterDelta() changes, no new sign convention.

Each outer iteration therefore becomes: (1) **exact bilateral step** --
solve :math:`S_{bb}\,\delta = r_{\text{old}}` via a precomputed factorization
and scatter the result, zeroing every bilateral row's residual in one shot;
(2) **contact-only sweep** -- run the usual ``gauss_seidel``/``jacobi``/
``colored`` sweep restricted to just the contact blocks (bilateral rows never
enter the sweep at all under ``true_schur``). ``chaotic`` is excluded (its
state lives in an atomic array, and combining that with an exact joint solve
raises questions not tackled here); setting ``condensed.true_schur=true``
with ``condensed.sweep_mode=chaotic`` is silently ignored.

**Factorization**: ``S_bb``'s sparsity, in row-space, is exactly the
shared-body adjacency graph restricted to bilateral rows -- a path graph for
every current example scene except ``hangbridge``, which has real branching
nodes (tripod apexes / deck planks, degree 4-6). ``BlockSparseLDLT``
(``src/misc/block_sparse_ldlt.{hpp,cpp}``) is a small, generic, reusable
block-sparse LDLT factorization over a graph of dense blocks (never a
dense/global matrix), used via
``CondensedAssembler::buildBilateralFactorization()``. Elimination order is
computed internally via a greedy minimum-degree heuristic, not left to
chance: an early version used natural (creation) order, which is exactly
correct but **hung** on ``hangbridge`` -- natural order is not fill-reducing
for a branching graph, and the fill-in cascade made factorization
impractically slow (confirmed by measurement, not just theory). Minimum-degree
recovers exactly the chain order (zero fill-in) on every path-graph scene, so
nothing is given up on the common case to fix the uncommon one; a synthetic
270-node stress test mirroring ``hangbridge``'s rope+plank topology factors
in under a millisecond with it.

The minimum-degree order is also **cached** across steps
(``CondensedAssembler`` keys the cache on the bilateral graph's structure --
dims + which pairs of blocks are coupled, not the numeric ``Gii``/
``complianceDiag`` values, which change every step regardless): every current
scene's bilateral topology is static for its lifetime (constraints aren't
created or destroyed at runtime), so after the first step this always hits,
skipping the O(n^2) symbolic pass entirely. ``factorWithOrder()`` still runs
the full *numeric* factorization on the current step's values either way, so
a cache hit only ever changes which (still fully valid) elimination order is
used, never the correctness of the factorization. Measured ~14% faster
``Condensed Setup`` on ``wilberforce`` (327.9µs to 282.6µs average, 2000
calls) -- a real but modest win, since ``Condensed Setup`` itself is a small
fraction of that scene's wall-clock: for the default one-VTK-frame-per-step
config, ``Output Write`` is ~90% of total time, unrelated to anything in this
solver. Worth knowing if chasing wall-clock further on a similar scene --
``output.interval_steps`` matters more there than anything below.

The cache's correctness under a runtime structural change (a constraint
actually added mid-simulation -- not exercised by any other example scene,
all of which have static bilateral topology) is verified by
``examples/scenes/dynamic_constraint``: three point masses, one spring from
the start, a second spring added via ``updateScene()`` partway through the
run. Confirmed (via ``debug.pj=true``) that ``nSprings`` actually grows from 1
to 2 mid-run, and that ``condensed`` with ``true_schur=true`` (cache
exercised, invalidates and recomputes at the structural change) produces a
``totalKE``/``posNormSum`` fingerprint identical to both PGS and ``condensed``
with ``true_schur=false`` -- not just an argument that the cache *should*
invalidate correctly, a scene that actually forces it to.

.. warning::
   **This is not a free "always enable it" win, and its benefit is much
   narrower than an earlier version of this note claimed -- measure before
   trusting it on a new scene/config.** An earlier benchmark reported a
   "~25x speedup, ~2000/step to ~2/step" finding for slinky; that was an
   artifact of reading only a truncated tail of a longer run's log, which
   happened to show just the *first* simulation step -- before ``slinky``'s
   coiled chain has collapsed into self-contact. Re-measured properly (full
   step-by-step logs, multiple independent runs, cross-checked against a
   from-scratch build of the pre-``true_schur`` commit to rule out an
   unrelated regression), the real picture is:

   - **Domino** (no springs): ``true_schur`` is a verified no-op --
     bit-identical ``totalKE``/position fingerprint with it on vs off, as
     required by construction (there is nothing to eliminate). This part of
     the original claim still holds.
   - **Slinky, step-by-step**: at ``gauss_seidel``, ``alpha=0.3``, no
     acceleration, tight tolerance (``1e-8``), the *first* step (zero
     contacts yet -- pure spring-chain relaxation) needs 2 outer iterations
     with ``true_schur=true`` vs 2532 without -- a real, large, reproducible
     effect confirming the theory for the case it actually applies to: no
     contacts at all. But ``slinky``'s coiled geometry means self-contacts
     start forming by the *second* step (414 contacts) and keep growing
     (500 → 800 → 1200 → 1500 → 2000+ within the first ~10 steps). From that
     point on, iteration counts are essentially the same with
     ``true_schur`` on or off (e.g. step 7: 16393 vs 16243; step 13: 46417
     vs 46328) -- contact active-set resolution, which ``true_schur`` does
     not touch, dominates almost immediately once contacts exist. Over a
     full 20-step run this scene spends only one step in the regime
     ``true_schur`` helps, so the *net* effect over the run was a **slowdown**
     (109.8s vs 85.4s) -- the per-iteration ``LDLT`` solve cost is paid on
     every step, but only repaid on the (rare, here) contact-free ones.
   - **Hangbridge** (branching topology, contact-heavy from the start -- 232
     contacts on 441 bodies): correctness verified (``totalKE`` agrees with
     PGS to within convergence-tolerance-level precision, ~1e-6 relative),
     but performance was roughly a wash to slightly worse -- consistent with
     the slinky finding above: this scene has contacts (and therefore the
     bottleneck ``true_schur`` doesn't address) from the very first step.

   **Practical guidance, corrected**: ``condensed.true_schur=true`` is a
   real, large win specifically for the portion of a simulation with **no
   active contacts** -- e.g. a compliant structure settling under gravity
   before it touches anything, or between contact events. It is not a
   general contact-scene speedup, and for a scene that is in contact
   (self- or otherwise) for most of its runtime -- which includes ``slinky``
   almost immediately -- expect it to be neutral-to-slower over a full run,
   not the dramatic win the first version of this section claimed. Measure
   iteration counts *per step*, not just a run's total or its first entry,
   before trusting a number on a new scene.

Nesterov acceleration
----------------------

``gauss_seidel``, ``jacobi``, and ``colored`` all support Nesterov (FISTA-style)
momentum acceleration via the standard ``pj.nesterov`` / ``pj.nesterov_beta_threshold``
/ ``pj.nesterov_restart_limit`` keys -- the exact same ones PJ and PGS already
read, reused rather than duplicated. It's implemented once, generically: a
``doSweep(lambda, u_corr)`` closure dispatches to whichever sweep mode is
configured, and the Nesterov loop only ever extrapolates/restarts the outer
``(lambda, u_corr)`` state between calls to it -- it has no knowledge of, and
requires no changes to, the sweep internals. ``chaotic`` is excluded: its
state lives in an atomic array rather than a plain vector, and
momentum-extrapolating a value that is *already* being deliberately kept
stale by design raises questions not addressed here.

.. warning::
   ``jacobi`` + Nesterov is not robust on every scene, but **not for the reason
   an earlier version of this note claimed**. It converges cleanly on
   ``domino`` (reaching the same fixed point as PJ's own Nesterov-accelerated
   run, in essentially the same sweep count), but diverges to `nan` on
   ``slinky`` at that scene's default ``pj.alpha=0.3``. This is **not** Nesterov
   shrinking Jacobi's stability margin -- a plain, unaccelerated ``jacobi``
   sweep (Nesterov off entirely) diverges on ``slinky`` at the same alpha, so
   momentum was never the cause. It is the same per-sweep-mode instability
   already documented in ``domino``'s ``scene_condensed.config``: ``jacobi``
   has a smaller stability margin than ``gauss_seidel``/``colored`` at a given
   ``pj.alpha``, and ``slinky``'s default alpha is tuned for the latter.
   Dropping ``pj.alpha`` to ``0.02`` makes ``jacobi`` + Nesterov stable on
   ``slinky`` too (confirmed over 1000+ steps with no divergence). The fix for
   a new scene's ``jacobi`` divergence is therefore to retune ``pj.alpha``
   first, not to disable Nesterov.

   Separately, the Nesterov loop's restart-on-non-finite-residual branch was
   at one point dead code: an early implementation ran every residual check
   through a throwing helper, so a `nan` fired an uncaught
   ``std::runtime_error`` before the restart logic ever got a chance to run,
   even though that logic exists precisely to recover from this case (ported
   from PJ's own ``nesterov_loop()``, which never throws). This has been
   fixed: the Nesterov loop now computes the residual without throwing, lets
   the existing restart/momentum-disable logic react to a non-finite value
   the same way it reacts to a merely-growing one, and only raises
   ``std::runtime_error`` (matching PJ/PGS's inherited behavior) once
   momentum has been fully disabled and the plain sweep output itself is
   still non-finite -- i.e. once there is genuinely nothing left to restart
   to. The extrapolation step also now checks the freshly-extrapolated state
   for finiteness before it is ever fed into the next sweep (``betak1`` is
   bounded to ``[0,1]``, but the extrapolated vector's magnitude is not, so
   this is a distinct failure mode from a bad ``betak1``).

Config keys
-----------

Iteration count, tolerance, relaxation, ``alpha``, warmstart, and Nesterov are
the same shared knobs PGS already reuses -- see :doc:`projected_jacobi`'s
config table for ``solver.max_iterations`` / ``solver.tol_abs`` /
``solver.tol_rel`` / ``solver.relaxation`` / ``solver.alpha`` /
``solver.warmstart`` / ``solver.nesterov``.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Key
     - Effect
   * - ``condensed.sweep_mode``
     - ``gauss_seidel`` (default) | ``jacobi`` | ``colored`` | ``chaotic``
   * - ``condensed.local_solve``
     - ``projection`` (default) | ``newton``
   * - ``condensed.num_threads``
     - OpenMP thread count for ``jacobi``/``colored``/``chaotic`` (0 = OpenMP default). More threads is not always faster -- see below.
   * - ``condensed.newton_max_iters``
     - Newton iteration cap per contact per call (default 8)
   * - ``condensed.newton_tol``
     - Newton convergence tolerance on :math:`\|\Phi\|` (default 1e-10)
   * - ``condensed.chaotic_reshuffle_interval``
     - Sweeps between block-order reshuffles in ``chaotic`` mode (default 50)
   * - ``condensed.chaotic_seed``
     - RNG seed for ``chaotic`` mode's shuffles
   * - ``condensed.true_schur``
     - Exactly eliminate bilateral (spring+damper) rows every outer iteration instead of sweeping them (default ``false``). No effect with ``sweep_mode=chaotic``. See `True Schur-complement elimination of the compliant chain (condensed.true_schur)`_ above -- validated as a large win in one regime and a wash/regression in another, not a blind default-on.

Known performance characteristics
-----------------------------------

Measured on ``examples/scenes/domino`` (pure rigid, ~13k frictional
contacts) and ``examples/scenes/slinky`` (compliant chain + friction, sparser
contact set); see ``CONDENSED_SOLVER_REPORT.md`` for the full numbers,
methodology, and caveats (short runs, one machine, not an exhaustive sweep).
These are specific, reproducible findings from that comparison, not general
claims:

- The per-block hot path (residual, update, scatter) uses fixed-size
  ``Vectorr<6>`` stack buffers, not dynamically-sized ``VectorXr`` -- an
  earlier version allocated several small ``VectorXr`` temporaries per block
  per sweep, which dominated wall-clock on scenes with many contacts. This
  closed roughly a third of the wall-clock gap to PGS on ``gauss_seidel``
  (domino: 19.16s to 13.64s for a 3-step run) and more on ``jacobi`` (8.37s
  to 5.18s), with zero change in converged answer or sweep count (verified
  bit-identical on the deterministic modes before/after).
- **With Nesterov enabled, condensed can beat PJ**: on domino, ``jacobi`` +
  Nesterov at 8 threads reaches ~1.8x PJ's speed (3.02s vs. PJ's 5.35s on a
  clean sequential run); on slinky, ``gauss_seidel`` + Nesterov reaches ~1.4x
  PJ's speed (8.39s vs. 12.04s). See the warning under `Nesterov
  acceleration`_ above -- ``jacobi`` + Nesterov is not safe on every scene.
- ``newton`` cuts iteration count roughly 10x on slinky at a single
  well-converged step. Without Nesterov, it remains faster in wall-clock over
  a 200-step run despite the extra per-block cost; layered *on top of*
  Nesterov (which already gets sweep count very low), the extra per-contact
  Newton cost was not repaid in the one comparison run so far.
- ``colored`` currently loses to ``jacobi``/``gauss_seidel`` on domino
  despite being parallel: ~50 colors x thousands of sweeps means ~50x the
  OpenMP barrier count of ``jacobi``'s 2-barriers-per-sweep, and that
  overhead dominates at this scale. Not recommended as a default until the
  coloring is cached across steps and/or more work is amortized per barrier.
- More threads is not always better, and the sweet spot can shift: on
  domino's ``jacobi`` mode a full thread-count sweep (1/2/4/8/16, sequential
  runs) found 8 threads fastest on this 16-core machine, with 16 threads
  measurably *slower* than 8 -- parallelization/scheduling overhead exceeding
  benefit at this per-thread work granularity. Tune ``condensed.num_threads``
  per scene rather than leaving it at the default of 0.
