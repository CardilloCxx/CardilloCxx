Projected Jacobi
================

The default solver. It handles both compliant constraints and non-smooth contact
simultaneously. The key idea is to decouple the contact impulse iterations from
the constraint solve using an LU-factorized free-motion system that is cheap to
apply repeatedly.

The problem
-----------

The solver minimises the kinetic energy deviation from the unconstrained
(free-motion) step subject to contact feasibility. It solves the scaled
system introduced in the Moreau chapter :doc:`../moreau_time_stepping`:

.. math::

   \mathcal{S}\,\mathbf{x}_{n+1} = \mathbf{b}_n,

where :math:`\mathcal{S}` is the scaled system matrix and :math:`\mathbf{b}_n` is the
assembled right-hand side (see the Moreau chapter for the exact block
structure, the position-level force law and the trajectory source-velocity
contributions). Contact impulses :math:`\Lambda_{NT}` enter the velocity block
of :math:`\mathbf{b}_n` via the term :math:`\mathbf{W}_{NT}\,\Lambda_{NT}`.

Fixed-point formulation and iteration
-------------------------------------

We treat contact impulses :math:`\mathbf{\Lambda}_{NT}` as the primary unknowns. The 
system is governed by the Delassus operator :math:`\mathbf{G}`:

.. math::

   \mathbf{G} = \mathbf{W}_{NT}^\top\,\mathcal{S}^{-1}\,\mathbf{W}_{NT}

This allows us to express the contact-space velocity :math:`\boldsymbol{\gamma}` as 
an affine map of the impulses:

.. math::

   \boldsymbol{\gamma}(\mathbf{\Lambda}_{NT}) = \mathbf{G}\,\mathbf{\Lambda}_{NT} + \mathbf{v}_{NT}

where :math:`\mathbf{v}_{NT}` is the contact-space velocity induced by 
movement of non-physics bodies and motors.

The Projected Jacobi iteration updates the impulses :math:`\mathbf{\Lambda}_{NT}` 
by moving toward the root of this map and projecting onto the feasible contact set 
:math:`\mathcal{P}`:

1. **Jacobi Sweep**:
   Update the current guess :math:`\mathbf{\Lambda}_{NT}^{(k)}` using the 
   diagonal step-size matrix :math:`\mathbf{R} \approx \text{diag}(\mathbf{G})^{-1}`:

   .. math::

      \mathbf{\Lambda}_{NT}^{(*)} = \mathbf{\Lambda}_{NT}^{(k)} - \mathbf{R}\,\boldsymbol{\gamma}(\mathbf{\Lambda}_{NT}^{(k)})

2. **Projection**:
   Project the result onto the feasible set :math:`\mathcal{P}` (Signorini for 
   normals, Coulomb friction for tangents):

   .. math::

      \mathbf{\Lambda}_{NT}^{(k+1)} = \mathcal{P}\left(\mathbf{\Lambda}_{NT}^{(*)}\right)

3. **Velocity Update**:
   Update the system state based on the change in impulses:

   .. math::

      \mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \mathcal{S}^{-1}\,\bar{\mathbf{W}}_{NT}\,\left(\mathbf{\Lambda}_{NT}^{(k+1)} - \mathbf{\Lambda}_{NT}^{(k)}\right)

The implementation warm-starts :math:`\mathbf{\Lambda}_{NT}^{(0)}` from the 
previous time step. It can be shown that the fixed point of this iteration 
satisfies the feasibility conditions for contact complementarity 
and Coulomb friction.


How it solves
-------------

.. mermaid::

   flowchart TD
       A[Assemble and factorise S with SparseLU] --> B[Compute free-motion velocity v_free]
       B --> C{Contacts?}
       C -- No --> D[Return v_free]
       C -- Yes --> E[Warmstart impulses from previous step]
       E --> F[Jacobi sweep: update all impulses in parallel]
       F --> G[Project impulses onto contact cone]
       G --> H[Convert impulse delta to velocity correction via S inverse]
       H --> I{Converged?}
       I -- No --> F
       I -- Yes --> J[Store impulses for warmstart, return v]

**Free-motion solve**
   S is assembled once per step and factorised with a sparse LU. Its sparsity
   pattern, and with it the symbolic analysis of the factorization
   (``analyzePattern``), is kept between steps: only the values change, until an
   entry becomes nonzero for the first time or bodies, constraints or used rows
   change (``misc::FixedPatternMatrix``, also used for the spring and damper rows
   :math:`W_g`, :math:`W_\gamma` of the ``DynamicsAssembler``). The free-motion
   velocity is the solution when no contacts are active. If there are no contacts
   the solver returns immediately.

**Jacobi sweep**
   Each iteration updates all impulse components simultaneously (in parallel,
   unlike Gauss-Seidel). The step size is scaled per row by the diagonal of the
   Delassus operator ``W M^-1 W^T``.

**Feasability projection**
   After each sweep, contact impulses are projected to satisfy:

   - **frictionless**: :math:`\Lambda_{N} \geq 0`
   - **frictional**: :math:`\|\Lambda_{T}\| \leq \mu\,\Lambda_{N}`

**Convergence**
   Checked every 10 iterations on the relative velocity residual
   (``velocity_rms``), or in every sweep on the natural map of the contact law
   (``solver.criterion = natural_map``), see *Convergence check* below. Controlled by
   ``solver.tol_abs``, ``solver.tol_rel``, and ``solver.max_iterations``.

Optional Nesterov acceleration
-------------------------------

When ``solver.nesterov = true``, momentum is added to the impulse updates using
the Nesterov scheme. The momentum coefficient is computed from the iteration
count and automatically reset when:

- the residual increases by more than 5%
- the momentum coefficient exceeds ``solver.nesterov_beta_threshold``
- the velocity update direction reverses

After ``solver.nesterov_restart_limit`` resets, momentum is disabled for the
remainder of the solve.

Warmstart
---------

When ``solver.warmstart = true``, the contact impulses from the previous step
are used as the initial guess. The contact tracker matches contacts across steps
by proximity (within ``collision.match_max_dist``) to identify which stored
impulse belongs to which current contact.


Convergence check
-----------------

We use a mixed absolute/relative residual:

.. math::

   \mathrm{res}^{(k)}
   = \frac{1}{\sqrt{n_u}}
     \left\lVert
       \frac{u^{(k+1)} - u^{(k)}}
            {\varepsilon_\mathrm{abs}\,\mathbf{1} + \varepsilon_\mathrm{rel}\,\max(|u^{(k+1)}|,|u^{(k)}|)}
     \right\rVert_2 \leq 1

This ``velocity_rms`` criterion is the default. It measures an increment, not
an error: averaged over all degrees of freedom, a single contact far from
convergence barely changes it, and with a small step ``solver.alpha`` the
increment is much smaller than the remaining error. Its tolerance therefore
does not bound the accuracy of the contacts.

With ``solver.criterion = natural_map``, PJ stops on the natural map of the contact
law with its own step :math:`\alpha R` (see :doc:`index`, section *Accuracy
measure*). The increment of the projection divided by the step,

.. math::

   r_b = (\alpha R_b)^{-1}\bigl(\Lambda_b^{(k)} - \Lambda_b^{(k+1)}\bigr),

is the natural map at the input of the sweep, a velocity error, and PJ stops
when :math:`e \le 1` with ``pj.tol_abs`` and ``pj.tol_rel`` in velocity units.
Divided by the full step, it does not shrink with ``solver.alpha`` as the
increment does: for a sticking contact it is the velocity error for any
:math:`\alpha`. It is computed inside the projection loop at a few operations
per contact, so it is evaluated in every sweep, also in the standard loop. The
final check (``solver.final_check = natural_map``) measures the result with the
common step of all solvers instead; near convergence the two agree.

Config keys
-----------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Key
     - Effect
   * - ``solver.max_iterations``
     - Hard iteration cap (default 200000)
   * - ``solver.tol_abs``, ``solver.tol_rel``
     - Tolerances of the selected criterion (``natural_map``: m/s)
   * - ``solver.criterion``
     - ``velocity_rms`` (default) or ``natural_map``: stop on the velocity increment or on the natural map of the contact law
   * - ``solver.final_check``
     - ``none`` (default) or ``natural_map``: natural map of the result after every step (all solvers)
   * - ``solver.alpha``
     - Step-size multiplier for each Jacobi sweep (default 0.3)
   * - ``solver.nesterov``
     - Enable Nesterov momentum (default false)
   * - ``solver.warmstart``
     - Reuse previous-step impulses as initial guess (default true)
