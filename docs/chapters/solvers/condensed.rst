Condensed (Block-Sparse) Solver
================================

.. contents:: On this page
   :local:
   :depth: 3

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

On top of that, ``condensed.true_schur`` can eliminate the compliant
(spring+damper) rows *exactly*, in closed form, instead of letting them
participate in the sweep at all -- see `True Schur-complement elimination of
the compliant chain (condensed.true_schur)`_ below.

Motivation: unlike :doc:`projected_jacobi` (PJ), which factorizes the full
mass+bilateral system once but re-solves it (a global sparse triangular
back-substitution) on **every** sweep, ``condensed`` -- like PGS -- never
touches a global system matrix after setup. Unlike PGS, it is not hardwired
to one sequential order or one local-solve rule, so the sweep/local-solve
choice can be tuned per scene, and the sweep loop can use OpenMP.

**Acceleration:** ``gauss_seidel``, ``jacobi``, and ``colored`` support both Nesterov
(``pj.nesterov``) and Chebyshev semi-iterative (``pj.chebyshev``) acceleration, the
same config keys PJ and PGS already use (no new keys) -- see `Nesterov
acceleration`_ and `Chebyshev acceleration`_ below. Anderson (PJ-only today) is
not yet ported.

From the Moreau system to condensed's row-blocks
--------------------------------------------------

This section derives *why* the ``RowBlock`` fields in the next section have
the values they have. It is not an independent re-derivation: every formula
below is a direct transcription of ``CondensedAssembler``'s actual C++ (cross-
checked line by line against ``PgsAssembler``, and validated empirically --
bit-identical results to PGS on every non-``true_schur``, non-gyroscopic
config). Where this section says a formula "is" something, that is a
statement about the code, not an independent mathematical claim.

**The starting point** is the scaled block linear system every Moreau-based
solver in this codebase solves, derived in full in
:doc:`../moreau_time_stepping`:

.. math::

   \mathcal{S}(\mathbf{q}_{n+\theta}, \Delta t, \theta)\;\mathbf{x}_{n+1} = \mathbf{b}_n

with :math:`\mathbf{x}_{n+1} = (\mathbf{u}_{n+1}, \boldsymbol\Lambda_{g,n+1}, \boldsymbol\Lambda_{\gamma,n+1})`
(contacts appended as extra rows/columns with a zero compliance block -- see that
chapter for the full block layout, and the note on :math:`\mathbf{M}_\mathrm{eff}`
below).

**Eliminating the velocity** (matching :doc:`projected_gauss_seidel`'s own
derivation, generalized to :math:`\mathbf{M}_\mathrm{eff}` instead of
:math:`\mathbf{M}`) gives a reduced system purely in the stacked constraint/contact
multipliers :math:`\boldsymbol\Lambda = (\boldsymbol\Lambda_g, \boldsymbol\Lambda_\gamma, \boldsymbol\Lambda_N, \boldsymbol\Lambda_T)`
(springs, then dampers, then frictionless contacts, then frictional contacts --
this exact packed order is ``CondensedTopology``'s own row layout):

.. math::

   \hat S\,\boldsymbol\Lambda = \tilde b,
   \qquad
   \hat S = \hat W^\top M_\mathrm{eff}^{-1} \hat W + \hat C,
   \qquad
   \hat W = (W_g, W_\gamma, W_N, W_T)

with :math:`\hat C` block-diagonal (nonzero only on the spring/damper rows --
contact rows have no compliance). Once :math:`\boldsymbol\Lambda` is known:

.. math::

   \mathbf{u}_{n+1} = \underbrace{M_\mathrm{eff}^{-1}\bigl(\mathbf{M}\mathbf{u}_n + \Delta t\,\mathbf f^\mathrm{ext}_{n+\theta}\bigr)}_{u_\mathrm{free}}
                    - \underbrace{M_\mathrm{eff}^{-1}\,\hat W\,\boldsymbol\Lambda}_{u_\mathrm{corr}}

``condensed`` -- like PGS -- never assembles :math:`\hat W` or :math:`\hat S` as
matrices (dense or sparse). Instead it decomposes :math:`\hat W`'s **row space**
into small, independent ``RowBlock`` s (one spring, one damper, one frictionless
contact, or one 3-row frictional contact each), and builds each block's own
slice of :math:`\hat S`/:math:`\tilde b` directly from two small dense Jacobians
and the bodies' (inverse) mass -- never touching a global matrix at any point.
The rest of this section gives the exact per-block formulas, in the same
notation ``condensed_assembler.cpp`` uses.

**Per-block local Delassus block** (``RowBlock::Gii``, built in
``CondensedAssembler::buildTopology()``) is exactly that block's own diagonal
block of :math:`\hat W^\top M_\mathrm{eff}^{-1}\hat W`:

.. math::

   G_i = J_a\,M_{\mathrm{eff},a}^{-1}\,J_a^\top + J_b\,M_{\mathrm{eff},b}^{-1}\,J_b^\top

:math:`J_a,J_b` are this block's own rows of :math:`\hat W` restricted to its
two touching bodies (``RowBlock::Ja``/``Jb``, zero-sized if that side is
static). :math:`M_{\mathrm{eff},a}^{-1}` is body :math:`a`'s own 3x3 (point
mass) or 6x6 (rigid body) block of :math:`M_\mathrm{eff}^{-1}` -- for the
overwhelming majority of bodies this is just ``diag(MinvDiag)`` (plain
per-DOF inverse mass), computed via Eigen's ``.asDiagonal()`` with no matrix
ever materialized; only a rigid body with an active implicit-gyroscopic
override (see `Implicit gyroscopic forces (moreau.implicit_gyroscopy)`_) gets a genuinely non-diagonal
6x6 block here, looked up from ``CondensedTopology::gyroMinvBlocks``.

**Per-block compliance** (``RowBlock::complianceDiag``, filled by
``CondensedAssembler::updateCompliance()``) is that block's own diagonal
slice of :math:`\hat C`: :math:`C_\mathrm{row}/(\theta\,\Delta t^2)` for
springs, :math:`A_\mathrm{row}/(\theta\,\Delta t)` for dampers, zero for
contacts -- read straight off ``DynamicsAssembler::Cdiag()``/``Adiag()``, no
assembly needed since :math:`\hat C` is block-diagonal by construction.

**Per-block right-hand side** (``CondensedAssembler::rhs()``, building
:math:`\tilde b` block by block) differs by kind, matching
:doc:`../moreau_time_stepping`'s :math:`\mathbf b_n` variable for variable
(:math:`\mathbf v_g \to` ``C_v_vec``, :math:`\mathbf v_\gamma \to`
``A_v_vec``, :math:`\mathbf g(t,\mathbf q) \to` ``g_error_vec``,
:math:`\beta \to` ``constraint_bias_factor``), in this codebase's own
internal sign convention (identical to ``PgsAssembler::rhs()`` -- not
independently re-derived here):

.. math::

   \text{rhs}_\mathrm{vel} = M\,\mathbf u_n + \Delta t\,\mathbf f^\mathrm{ext} \;\; [+\, \Delta t\,\mathbf f^\mathrm{gyr} \text{ if gyroscopic forces are NOT treated implicitly}]

For a **spring** row (:math:`\Lambda_{g}` held at its warm-started/previous-step
value :math:`\Lambda_{g,n}` while solving for the current step's answer --
this is a fixed constant *within* one ``solve()`` call, not the unknown being
solved for):

.. math::

   \mathrm{seg}_\mathrm{spring} = \frac{1}{\theta\Delta t^2}\,C_\mathrm{row}\odot\Lambda_{g,n}
                                 + \frac{1-\theta}{\theta}\bigl(J_a\mathbf u_{n,a}+J_b\mathbf u_{n,b}\bigr)
                                 + \frac{1}{\theta}\,\mathbf v_{g,\mathrm{row}}
                                 + \bigl[\beta>0\bigr]\,\gamma_\mathrm{stab}
                                 + J_a M_{\mathrm{eff},a}^{-1}\,\text{rhs}_{\mathrm{vel},a} + J_b M_{\mathrm{eff},b}^{-1}\,\text{rhs}_{\mathrm{vel},b}

with :math:`\gamma_\mathrm{stab} = \bigl(-C_\mathrm{row}\odot\Lambda_{g,n}/\Delta t + \mathbf g_\mathrm{row}\bigr)\,\beta/(\Delta t\,\theta)`.
A **damper** row drops the compliance/bias terms (dampers carry no position
error to stabilize):

.. math::

   \mathrm{seg}_\mathrm{damper} = \frac{1-\theta}{\theta}\bigl(J_a\mathbf u_{n,a}+J_b\mathbf u_{n,b}\bigr)
                                 + \frac{1}{\theta}\,\mathbf v_{\gamma,\mathrm{row}}
                                 + J_a M_{\mathrm{eff},a}^{-1}\,\text{rhs}_{\mathrm{vel},a} + J_b M_{\mathrm{eff},b}^{-1}\,\text{rhs}_{\mathrm{vel},b}

A **contact** row (frictionless or frictional) is simplest of all -- it only
ever needs :math:`u_\mathrm{free}`, already fully formed, plus the
kinematic contact-velocity bias:

.. math::

   \mathrm{seg}_\mathrm{contact} = J_a\,u_{\mathrm{free},a} + J_b\,u_{\mathrm{free},b} + \mathbf b_{\mathrm{contact},\mathrm{row}}

.. note::
   Springs/dampers and contacts use two **formally different, not merely
   differently-named** quantities for the "free velocity" term:
   :math:`J\,M_\mathrm{eff}^{-1}\,\text{rhs}_\mathrm{vel}` for bilateral rows,
   :math:`J\,u_\mathrm{free}` for contact rows. They coincide only when
   ``moreau.lambda_theta=false`` (the default): with it enabled, an extra
   correction term :math:`-(1-\theta)\,\hat W^\top(\ldots\Lambda_{g,n}\ldots\Lambda_{\gamma,n}\ldots)`
   is folded into :math:`\text{rhs}_\mathrm{vel}` (an alternate theta-scaling
   of the constraint-force contribution) but **not** into
   :math:`u_\mathrm{free}` (built independently by ``CondensedAssembler::ufree()``,
   which has no knowledge of ``moreau_lambda_theta`` at all) -- so bilateral
   and contact rows genuinely diverge in that mode, by design, not by
   oversight.

Building :math:`\hat S`'s off-diagonal coupling (needed only by
``condensed.true_schur``, since every sweep mode above only ever needs each
block's own :math:`G_i`, never a cross-block term) is described in
`How Sbb's off-diagonal blocks are built`_ below.

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
     - :math:`G_i = J_a\,M_{\mathrm{eff},a}^{-1}\,J_a^\top + J_b\,M_{\mathrm{eff},b}^{-1}\,J_b^\top` -- see `From the Moreau system to condensed's row-blocks`_
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

.. note::
   :math:`M_{\mathrm{eff}}^{-1}` above is *usually* just ``diag(MinvDiag)`` --
   a plain per-DOF inverse mass, applied via Eigen's ``.asDiagonal()``
   (a column-scaling operation, no matrix ever materialized). Every one of
   ``Gii``'s two accumulation terms, and the analogous terms in ``rhs()``/
   ``ufree()``/the bilateral coupling below, individually check whether the
   relevant body has an active implicit-gyroscopic override
   (``CondensedTopology::gyroMinvBlocks``) and only pay the cost of a real
   6x6 block-times-block product on that (rare) path -- see `Implicit
   gyroscopic forces (moreau.implicit_gyroscopy)`_.

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
:math:`r_i` the block's local residual, :math:`D_i^{-1}` its ``GiiInv``) --
they differ only in *when* each block reads/writes the shared
velocity-correction state ``u_corr``.

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
   sweep. See `Graph coloring algorithm`_ below for exactly how the coloring
   is computed.

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

Graph coloring algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~

``colorGreedyWelshPowell()`` (``src/misc/graph_coloring.hpp``/``.cpp``) is a
standard greedy graph coloring, computed fresh every ``solve()`` call when
``condensed.sweep_mode=colored`` (unlike ``true_schur``'s elimination order,
this is **not** cached across steps today -- see the performance note below):

1. **Adjacency**: two blocks are adjacent iff they share a body
   (``CondensedTopology::blocksOfBody``, restricted to the contact-block
   range when ``condensed.true_schur`` is active -- bilateral rows never
   enter the sweep in that case, so they're excluded from this graph too).
2. **Order nodes by descending degree** (Welsh-Powell's own heuristic: color
   the most-constrained nodes first, when the fewest choices remain).
3. **Greedily assign each node the smallest color not used by any
   already-colored neighbor**, processed in that order. Implemented with a
   monotonically increasing "stamp" trick rather than a fresh
   ``O(numColors)`` reset per node: a ``forbiddenStamp[c] == stamp`` check
   marks color ``c`` as taken by the node currently being colored, and the
   stamp increments once per node -- so checking "is color c forbidden"
   is :math:`O(1)` regardless of how many colors exist so far.

The result (``Coloring::colorClasses``) is one contiguous block-index list
per color, consumed directly by ``coloredSweep()``: color 0's blocks are all
updated (in parallel) before any of color 1's blocks start, guaranteeing no
two concurrently-running blocks ever touch the same body's ``u_corr`` slice.

This produces more colors than a graph-theoretically optimal (chromatic-
number) coloring in general -- greedy Welsh-Powell is a well-known
polynomial-time heuristic, not an exact minimum-coloring algorithm, which is
NP-hard in general -- but is fast to compute (:math:`O(V\log V + E)`, the
sort dominating) and good enough in practice: see the performance note in
`Known performance characteristics`_ on why ``colored`` currently loses to
``jacobi``/``gauss_seidel`` on ``domino`` despite being parallel (too many
colors relative to the per-color OpenMP barrier cost at that scene's
contact density) -- caching this coloring across steps (like ``true_schur``'s
elimination order already is) is the natural next step if that overhead
ever becomes the bottleneck rather than the tuned-per-scene ``alpha``/sweep
count it currently is.

Local solve: projection vs. semismooth Newton
-----------------------------------------------

**projection** (default) clips exactly as :doc:`projected_gauss_seidel`
does: normal impulse clamped to :math:`\Lambda_N \le 0`, tangential impulses
scaled toward the Coulomb disk of radius :math:`-\mu\,\Lambda_N`.

**newton** replaces that clip, for frictional (3-row) blocks only, with a
guarded Alart-Curnier semismooth Newton solve
(``src/physics/solver/local_contact_newton.{hpp,cpp}``). Given the block's
local Delassus matrix :math:`G_i` and current residual :math:`r` (evaluated
at the current ``lambda``, i.e. ``blockResidual()``'s return value for this
block):

.. math::

   \rho_N = 1/G_{i,00},
   \qquad
   \rho_T = 1/\lambda_\mathrm{max}(G_{i,TT})

(:math:`G_{i,TT}` the 2x2 tangential sub-block; :math:`\lambda_\mathrm{max}`
via the closed-form 2x2 eigenvalue
:math:`\tfrac12\bigl((a+d)+\sqrt{(a-d)^2+4b^2}\bigr)`, no
``SelfAdjointEigenSolver`` needed -- :math:`\rho_N`/:math:`\rho_T` are
diagonal preconditioners for the local Newton iteration, analogous in role
to ``GiiInv`` for the linear path, just derived to make the *Newton* map
well-scaled rather than to directly precondition a linear solve). With
:math:`M = I - \mathrm{diag}(\rho_N,\rho_T,\rho_T)\,G_i` and
:math:`c = \mathrm{diag}(\rho)\,(r + G_i\,\lambda)`, the trial state
:math:`y(\lambda') = M\lambda' + c` (built so that :math:`y(\lambda_\mathrm{cur})`
reduces to exactly the plain relaxation update's trial state before
projection) is case-split, per Newton iterate, into:

- **separation** (:math:`y_N \ge 0`): no contact force, :math:`g=0`,
  :math:`\partial g/\partial\lambda = 0`.
- **sticking** (:math:`\|y_T\| \le -\mu y_N`): the unconstrained trial state
  is already inside (or on) the friction cone, :math:`g=y`,
  :math:`\partial g/\partial\lambda = M`.
- **sliding** (otherwise): the tangential trial state is projected onto the
  disk of radius :math:`-\mu y_N`, :math:`g_T = -\mu y_N\,\hat t`
  (:math:`\hat t = y_T/\|y_T\|`), with the Jacobian's tangential rows
  accounting for both the direct scaling and the trial state's own
  dependence on :math:`\lambda` through :math:`y_N` and :math:`\hat t`.

This builds :math:`\Phi(\lambda) = \lambda - g(\lambda)` and its Jacobian
:math:`J_\Phi = I - \partial g/\partial\lambda`, then takes a Newton step
:math:`J_\Phi\,\delta = -\Phi` (solved via ``Eigen::FullPivLU`` on the 3x3
:math:`J_\Phi`, checked for invertibility every iteration since the active
case -- and therefore the Jacobian's structure -- can change step to step).

The solve is **guarded** at every point that could otherwise silently
corrupt ``lambda`` -- it reports failure (leaving ``lambda`` completely
unchanged) on any of:

- :math:`G_{i,00} \le \varepsilon` or :math:`\lambda_\mathrm{max}(G_{i,TT}) \le \varepsilon`
  (near-singular local block -- :math:`\rho_N`/:math:`\rho_T` would blow up);
- :math:`J_\Phi` not invertible (``FullPivLU::isInvertible()``);
- a Newton step :math:`\delta` that is not finite;
- a step that *increases* :math:`\|\Phi\|` relative to the previous
  iterate, once past the very first iteration (the first step is always
  accepted even if it doesn't immediately reduce the residual, since the
  case-split can legitimately jump on iteration 0 before settling);
- exhausting ``condensed.newton_max_iters`` without meeting
  ``condensed.newton_tol`` on :math:`\|\Phi\|`.

On failure the caller falls straight through to the ``projection`` update
for that block -- ``newton`` never risks being *less* robust than
``projection``, only sometimes not worth its extra cost.

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

**The math**, worked from the same :math:`\hat S = \hat W^\top M_\mathrm{eff}^{-1}\hat W + \hat C`
system derived above (nothing below is re-derived independently of that
verified ground truth). Partition rows into bilateral (springs +
dampers, pure equality constraints -- ``projectBlock()`` already no-ops for
these) and contact (subject to the friction-cone projection):

.. math::

   \hat S = \begin{bmatrix} S_{bb} & S_{bc} \\ S_{cb} & S_{cc} \end{bmatrix},
   \quad
   \Lambda = \begin{bmatrix} \Lambda_b \\ \Lambda_c \end{bmatrix}

Eliminating :math:`\Lambda_b` gives the Schur complement over contacts alone:
:math:`(S_{cc} - S_{cb} S_{bb}^{-1} S_{bc})\,\Lambda_c = \text{rhs}_c - S_{cb} S_{bb}^{-1}\,\text{rhs}_b`.
Explicitly forming this (generally dense, since :math:`S_{bb}^{-1}` is dense
even though :math:`S_{bb}` is sparse) matrix is exactly what this solver
family exists to avoid, and turns out not to be necessary: because bilateral
rows are pure equalities with no projection, whatever the current contact
impulses are, there is a single exact :math:`\Lambda_b` that zeroes every
bilateral row's residual simultaneously --
:math:`\delta\Lambda_b = S_{bb}^{-1}\,r_{\text{old}}`, where :math:`r_{\text{old}}`
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

.. note::
   Because a direct linear solve of :math:`S_{bb}` always lands on the same
   unique solution regardless of the state it started from, one exact
   bilateral step is mathematically sufficient to zero the bilateral
   residual for good -- a *second* call, given the first call's exact
   output, is provably a no-op. Any outer-iteration count above 1-2 you
   observe on a scene with no contacts is overwhelmingly the surrounding
   convergence-check/Nesterov bookkeeping confirming that nothing changed
   anymore, not the Schur solve itself needing repeated work -- see
   ``CONDENSED_SOLVER_REPORT.md`` for a worked example (``wilberforce``) with
   exact iteration-count accounting.

How Sbb's off-diagonal blocks are built
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`S_{bb}`'s diagonal blocks are already sitting in ``RowBlock::Gii``/
``complianceDiag`` -- no extra work needed (``buildBilateralFactorization()``
just adds them: ``diagBlocks[i] = Gii_i + diag(complianceDiag_i)``). The
off-diagonal *coupling* between two bilateral blocks :math:`p,q` exists only
if they share a body :math:`b` (an end of one spring/damper is the same
rigid body as an end of another), and is built the same way ``Gii`` itself
is -- directly from the two blocks' own Jacobian slices at the shared body,
never through a global matrix:

.. math::

   S_{bb}[p,q] \mathrel{+}= J_p\, M_{\mathrm{eff},b}^{-1}\, J_q^\top

summed over every body two bilateral blocks share (a pair can share *two*
bodies at once -- e.g. two independent springs between the same two rigid
bodies -- so contributions are accumulated into a map keyed by the
unordered block-index pair before being handed to ``BlockSparseLDLT``, which
expects at most one edge per pair). In the overwhelming common case
(:math:`M_{\mathrm{eff},b}^{-1}` diagonal), :math:`S_{bb}[p,q]` and
:math:`S_{bb}[q,p]` are exact transposes of each other, so only one
direction is computed and stored per pair -- ``BlockSparseLDLT`` derives the
other via transpose (its ``symmetric=true`` mode, the default). When body
:math:`b` has an active implicit-gyroscopic override, that is no longer
true (:math:`(J_p M^{-1} J_q^\top)^\top = J_q (M^{-1})^\top J_p^\top \ne
J_q M^{-1} J_p^\top` for a non-symmetric :math:`M^{-1}`), so **both**
directions are computed independently for the *entire* bilateral graph once
any pair in it is affected (an edge given in only one direction in
``BlockSparseLDLT``'s ``symmetric=false`` mode is treated as exactly zero in
the other, not "derive it") -- see `Implicit gyroscopic forces (moreau.implicit_gyroscopy)`_.

Block-sparse factorization algorithm (``BlockSparseLDLT``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``S_bb``'s sparsity, in row-space, is exactly the shared-body adjacency graph
restricted to bilateral rows -- a path graph for every current example scene
except ``hangbridge``, which has real branching nodes (tripod apexes / deck
planks, degree 4-6). ``BlockSparseLDLT`` (``src/misc/block_sparse_ldlt.{hpp,cpp}``)
is a small, generic, reusable block-sparse factorization over a graph of
dense blocks (never a dense/global matrix), used via
``CondensedAssembler::buildBilateralFactorization()``.

**Elimination order.** Computed internally via a greedy minimum-degree
heuristic, not left to chance: at each step, eliminate whichever
still-active node currently has the fewest active neighbors, then simulate
the resulting fill-in (every pair of that node's still-active neighbors
becomes connected, mirroring the real Schur update below) purely on
adjacency, before ever touching a matrix block. An early version used
natural (creation) order, which is exactly correct but **hung** on
``hangbridge`` -- natural order is not fill-reducing for a branching graph,
and the fill-in cascade made factorization impractically slow (confirmed by
measurement, not just theory). Minimum-degree recovers exactly the chain
order (zero fill-in) on every path-graph scene, so nothing is given up on
the common case to fix the uncommon one; a synthetic 270-node stress test
mirroring ``hangbridge``'s rope+plank topology factors in under a
millisecond with it. This order is also **cached** across steps (see the
performance discussion further below).

**Symmetric mode (block LDLT, the default -- used whenever no implicit-
gyroscopic body touches the bilateral graph).** For each pivot :math:`k` in
elimination order, with :math:`D_k` its current (possibly already-updated by
earlier eliminations) diagonal block:

.. math::

   D_k^{-1} \leftarrow \texttt{invertSmallSpd}(D_k)

For every still-active neighbor :math:`p` of :math:`k`:

.. math::

   L_{p,k} = \mathrm{block}(p,k)\,D_k^{-1}

Then, for **every pair** :math:`(p,q)` of :math:`k`'s still-active neighbors
(including :math:`p=q`, the diagonal update) -- this is the block
generalization of one elimination step of scalar sparse Cholesky:

.. math::

   \mathrm{block}(p,q) \mathrel{-}= L_{p,k}\,\mathrm{block}(k,q)

For :math:`p=q` this updates :math:`p`'s own diagonal block directly
(re-symmetrized after each update, :math:`D_p \leftarrow \tfrac12(D_p+D_p^\top)`,
to keep it numerically symmetric against floating-point drift); for
:math:`p\ne q` this either updates an existing coupling block or creates a
*fresh fill-in edge* if :math:`p,q` weren't already connected, storing the
transpose at :math:`(q,p)` (only one direction is ever computed -- the
other is always the transpose, since :math:`D_k` is symmetric here). Every
pivot's :math:`D_k^{-1}` and its list of :math:`L_{p,k}` blocks are kept
(``PivotFactor``), in elimination order, forming the stored factorization.

**Solving** :math:`S_{bb}\,x = \mathrm{rhs}` then proceeds in three passes
over the stored pivots:

1. **Forward substitution** (:math:`Ly=\mathrm{rhs}`, :math:`L` unit lower-
   triangular in elimination order): visiting pivots in elimination order,
   :math:`y_p \mathrel{-}= L_{p,k}\,y_k` for every :math:`p` in pivot
   :math:`k`'s stored neighbor list.
2. **Apply** :math:`D^{-1}`: :math:`x_k = D_k^{-1}\,y_k` for every pivot.
3. **Backward substitution** (:math:`L^\top x = D^{-1}y`, reverse
   elimination order): :math:`x_k \mathrel{-}= L_{p,k}^\top\,x_p` for every
   :math:`p` in :math:`k`'s stored neighbor list.

**Non-symmetric mode (block LU -- only when at least one bilateral pair's
coupling is genuinely non-symmetric, see `Implicit gyroscopic forces (moreau.implicit_gyroscopy)`_).**
:math:`\mathrm{block}(k,q)` and :math:`\mathrm{block}(q,k)` are no longer
transposes of each other, so nothing can be derived from the other and
**both** an :math:`L` list (rows below the pivot) and a separate :math:`U`
list (columns to the right, stored raw, no :math:`D^{-1}` applied) are kept
per pivot:

.. math::

   D_k^{-1} \leftarrow \texttt{invertSmallGeneral}(D_k),
   \qquad
   L_{p,k} = \mathrm{block}(p,k)\,D_k^{-1},
   \qquad
   U_{k,q} = \mathrm{block}(k,q)

The Schur update now runs over every **ordered** pair :math:`(p,q)` of
active neighbors (both :math:`(p,q)` and :math:`(q,p)` computed
independently, not one derived from the other via transpose):

.. math::

   \mathrm{block}(p,q) \mathrel{-}= L_{p,k}\,U_{k,q}

roughly double the symmetric path's Schur-update work, since nothing can be
skipped via symmetry. Solving mirrors this: forward substitution is
identical in structure to the symmetric case (:math:`L` is unit lower-
triangular either way), but backward substitution solves :math:`Ux=y`
directly with the stored raw :math:`U` blocks, applying :math:`D^{-1}`
**after** subtracting the :math:`U`-contributions rather than before:

.. math::

   x_k = D_k^{-1}\Bigl(y_k - \sum_q U_{k,q}\,x_q\Bigr)

This is exactly the reason the two modes are two separate code paths rather
than one branchy general algorithm: paying the non-symmetric cost
unconditionally -- on every scene, including the overwhelming majority that
never touch implicit gyroscopic forces -- would be a pure regression for no
benefit. ``BlockSparseLDLT`` picks the mode once, at ``build()`` time
(``symmetric=true`` default), based on whether ``buildBilateralFactorization()``
ever needed a directed edge (see the section above).

**Caching the elimination order.** ``CondensedAssembler`` keys a cache on
the bilateral graph's *structure* (dims + which pairs of blocks are coupled,
not the numeric ``Gii``/``complianceDiag`` values, which change every step
regardless): every current scene's bilateral topology is static for its
lifetime (constraints aren't created or destroyed at runtime), so after the
first step this always hits, skipping the :math:`O(n^2)` minimum-degree
symbolic pass entirely. ``factorWithOrder()`` still runs the full *numeric*
factorization on the current step's values either way, so a cache hit only
ever changes which (still fully valid) elimination order is used, never the
correctness of the factorization. Measured ~14% faster ``Condensed Setup``
on ``wilberforce`` (327.9µs to 282.6µs average, 2000 calls) -- a real but
modest win, since ``Condensed Setup`` itself is a small fraction of that
scene's wall-clock: for the default one-VTK-frame-per-step config, ``Output
Write`` is ~90% of total time, unrelated to anything in this solver. Worth
knowing if chasing wall-clock further on a similar scene --
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

The other half of the question -- does the cache actually *hit* when the
structure doesn't change, not just invalidate correctly when it does -- is
directly measurable via ``CARDILLO_DEBUG_SCHUR_CACHE=1`` (prints a HIT/MISS
line with running counts every call). Confirmed on ``wilberforce``: 999 hits
out of 1000 steps, missing only the first call. On ``dynamic_constraint``:
misses exactly at step 1 and at the step the second spring is added, hits
every other call (98/100). Both match the theoretical expectation exactly --
this isn't inferred from a wall-clock delta, it's counted directly.

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

Implicit gyroscopic forces (``moreau.implicit_gyroscopy``)
------------------------------------------------------------

``moreau.implicit_gyroscopy=true`` (a second-order, Kahan-consistent treatment of the Euler
equations for a rotating rigid body) is now genuinely supported by ``condensed`` -- it was silently
dropped before (the assembler printed a warning and applied nothing, giving wrong dynamics for any
rigid body with the flag set, not just a missing optimization).

The correction makes a rigid body's effective mass non-symmetric: ``M_eff = M - dt*Grot`` where
``Grot = 0.5*([I*omega]_x - [omega]_x*I)`` (the same formula ``projected_jacobi`` already used, and
the same :math:`\mathbf M_\mathrm{eff}` that appears throughout `From the Moreau system to
condensed's row-blocks`_ above). ``condensed``'s Schur-complement architecture needs this as an
explicit per-body inverse-mass block rather than a system-matrix correction, so a body with the flag
active (and actual rotational dof) gets a ``Minv`` override -- unchanged diagonal for translation, a
general (not necessarily symmetric) inverse for rotation -- consulted everywhere the solver would
otherwise use the plain diagonal ``MinvDiag``: every ``Gii`` accumulation in `Block-sparse
representation`_, ``ufree()``/``rhs()``'s free-velocity terms, and the bilateral coupling in `How
Sbb's off-diagonal blocks are built`_. Every one of those call sites falls through to the *original*
expression, unchanged, for every body without an active override, so a scene that never enables this
flag is byte-for-byte unaffected.

When the compliant (spring/damper) chain includes such a body, ``condensed.true_schur``'s
block-sparse factorization also switches from block-LDLT to the block-LU generalization described in
`Block-sparse factorization algorithm (BlockSparseLDLT)`_ -- roughly double the per-pivot cost, but
only for that graph, and only when a gyroscopic body actually participates in it; a scene with
gyroscopic bodies that never touch the compliant chain keeps the cheaper symmetric path.

.. note::
   Validated against ``projected_jacobi``, the trusted independent reference for this formula, on
   ``wilberforce`` (a torsional-vertical mode-coupled pendulum discretized as 360 coupled rigid
   segments -- the standard textbook exercise for this exact effect): with
   ``condensed.true_schur=true`` and ``condensed.sweep_mode=colored``, kinetic energy after 1000
   steps agreed with PJ to ~10 significant figures, and the position fingerprint matched exactly to
   every printed digit. Building that PJ reference surfaced an independent, pre-existing bug in
   ``projected_jacobi`` itself -- its solver never passed ``moreau_implicit_gyroscopy``/
   ``moreau_lambda_theta`` through to the system-matrix assembly, so PJ's own implicit treatment was
   also inactive until fixed. See ``CONDENSED_SOLVER_REPORT.md`` for the full validation writeup.

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

.. note::
   The reported iteration count (``iters=`` in the progress bar/logs, backed
   by ``SolverBase::lastIterations()``) counts exactly one call to
   ``doSweep()``/``exactBilateralStep()`` per unit -- this was previously
   inflated by one for the Nesterov branch specifically (``iter+2`` instead
   of ``iter+1``, inherited verbatim from ``ProjectedGaussSeidelSolver``'s
   own Nesterov branch, which had the identical off-by-one). Neither branch
   performs a second (e.g. forward+backward) sweep to justify a higher
   count -- checked directly: every loop variant in this codebase (PGS's own
   plain/Nesterov branches, PJ's ``nesterov_loop``/``chebyshev_loop``/
   ``anderson_loop``, ``condensed``'s own plain/chaotic branches) calls its
   sweep function exactly once per outer iteration. Fixed in both PGS and
   ``condensed``; purely cosmetic (only the progress-bar/log display was
   affected, never any convergence or control-flow decision), but relevant
   if you're using iteration counts to characterize how much work an exact
   ``true_schur`` solve actually needs (see the note in `True
   Schur-complement elimination of the compliant chain (condensed.true_schur)`_ above).

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

Chebyshev acceleration
------------------------

``gauss_seidel``, ``jacobi``, and ``colored`` also support Chebyshev semi-iterative
acceleration via ``pj.chebyshev`` -- the same key ``projected_jacobi`` already uses
(no condensed-specific key), with the same precedence rule: ``pj.nesterov`` takes
priority if both are set (``useChebyshev = pj_chebyshev && !pj_nesterov``). Unlike
Nesterov's residual-dependent momentum/restart logic, Chebyshev uses a **fixed**
extrapolation schedule
(:math:`\omega_1 = 1/(1-\rho^2/2)`, :math:`\omega_{k+1} = 1/(1-\tfrac14\rho^2\omega_k)`,
:math:`x_{k+1} = \omega_{k+1}(G(x_k) - x_{k-1}) + x_{k-1}`) driven entirely by one
number, :math:`\rho`: an estimate of the spectral radius of the *linearized* sweep
operator :math:`G` (i.e. with ``projectBlock()``/Newton disabled -- exactly the
nonlinearities that make a rigorous spectral-radius argument inapplicable in the
first place, the same caveat ``projected_jacobi``'s own Chebyshev implementation
carries).

**Estimating** :math:`\rho` (``estimateSpectralRadius()`` in ``condensed_solver.cpp``,
computed once per ``solve()`` call, before the main iteration loop) has to work
without ever assembling a matrix, unlike PJ's own estimator (which power-iterates
through ``PjAssembler::solveS()``, a real sparse solve). The key observation: with
the sweep's rhs forced to zero and projection/Newton disabled, one sweep step
becomes an *exact* linear map :math:`G(d) = A\,d` -- no affine offset, since the
rhs-dependent constant vanishes -- so plain power iteration (8 rounds, a
deterministic all-ones probe vector, matching PJ's own choice of a fixed rather
than randomized start) on this zero-rhs, linear-only variant of the *actual
configured sweep* (``gaussSeidelSweep``/``jacobiSweep``/``coloredSweep``, plus
``exactBilateralStep`` under ``condensed.true_schur`` -- reusing the real sweep
functions with a ``linearOnly`` flag threaded through, not a separate
reimplementation) converges to :math:`|A|`'s dominant eigenvalue magnitude, exactly
analogous to PJ's power iteration on its own (matrix-based) linearized operator.
Clamped to :math:`\le 0.995`, same reasoning as PJ: the :math:`\omega` recurrence
divides by :math:`1-\tfrac14\rho^2\omega`, only well-defined for :math:`\rho<1`.

.. note::
   Verified: on ``wilberforce`` (pure DAE, ``condensed.true_schur=true`` exactly
   eliminates every row already, so there's nothing to accelerate) Chebyshev
   reproduces the un-accelerated fixed point exactly, bit for bit. On ``domino``
   and ``slinky`` (frictional-contact-dominated, ``colored`` sweep), Chebyshev
   needed **fewer total sweeps than Nesterov on both** in the one comparison run
   so far (domino: 1415 vs. 1514 Nesterov vs. 1774 unaccelerated, ~20% fewer than
   no acceleration; slinky: 1403 vs. 1418 Nesterov vs. 1754 unaccelerated) --
   consistent with the "matches or beats Nesterov" finding already established for
   PJ, now reproduced independently on ``condensed``'s different (matrix-free,
   block-sweep) architecture. Not an exhaustive sweep (one machine, short runs,
   two scenes) -- measure before trusting these specific numbers on a new scene,
   same discipline as everything else in this file.

.. warning::
   Chebyshev inherits the exact same ``jacobi``-at-too-high-``alpha`` instability
   documented in the Nesterov warning above -- confirmed directly: ``jacobi`` +
   Chebyshev diverges on ``domino`` at ``pj.alpha=1.0`` (tuned for
   ``colored``/``gauss_seidel``), and so do plain ``jacobi`` and ``jacobi`` +
   Nesterov at that same alpha on that same scene -- this is ``jacobi``'s own
   stability margin, not something any acceleration method is expected to fix.
   Retune ``pj.alpha`` down first (as the Nesterov warning already advises) if
   ``jacobi`` diverges on a new scene, regardless of which acceleration (if any)
   is enabled.

   Getting this case to fail *cleanly* took one extra guard beyond porting PJ's
   ``chebyshev_loop()`` directly: PJ/Nesterov both check whether the **raw sweep
   output itself** (before any extrapolation) has gone non-finite, and throw
   immediately if so -- there's nothing left to fall back to or extrapolate from
   at that point, only an already-diverged sweep. An initial port of Chebyshev
   only guarded the *extrapolation's own* output (matching the `nan`-safety
   already needed for the linearization-breakdown case), which is not the same
   check -- without it, once the raw sweep itself started producing `nan` (e.g.
   from the ``jacobi``-instability case above), the fallback path kept feeding
   that same non-finite state back into the sweep every remaining iteration,
   spinning silently for the rest of ``pj.max_iterations`` instead of failing
   fast. Fixed by adding the identical raw-output check the Nesterov branch
   already had.

Config keys
-----------

Iteration count, tolerance, relaxation, ``alpha``, warmstart, Nesterov, and
Chebyshev are the same shared knobs PGS/PJ already reuse -- see
:doc:`projected_jacobi`'s config table for ``solver.max_iterations`` /
``solver.tol_abs`` / ``solver.tol_rel`` / ``solver.relaxation`` /
``solver.alpha`` / ``solver.warmstart`` / ``solver.nesterov`` / ``pj.chebyshev``.

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
  ``Vector6r`` stack buffers, not dynamically-sized ``VectorXr`` -- an
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
  coloring is cached across steps (see `Graph coloring algorithm`_ above)
  and/or more work is amortized per barrier.
- More threads is not always better, and the sweet spot can shift: on
  domino's ``jacobi`` mode a full thread-count sweep (1/2/4/8/16, sequential
  runs) found 8 threads fastest on this 16-core machine, with 16 threads
  measurably *slower* than 8 -- parallelization/scheduling overhead exceeding
  benefit at this per-thread work granularity. Tune ``condensed.num_threads``
  per scene rather than leaving it at the default of 0.
