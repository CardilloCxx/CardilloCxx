Projected Jacobi
================

The default solver. It handles both compliant constraints and non-smooth contact
simultaneously. The key idea is to decouple the contact impulse iterations from
the constraint solve using an LU-factorized free-motion system that is cheap to
apply repeatedly.

The problem
-----------

The solver minimises the kinetic energy deviation from the unconstrained
(free-motion) step subject to contact feasibility. Concretely it solves:

.. math::

   S \begin{pmatrix} u \\ \lambda_g \\ \lambda_\gamma \end{pmatrix} = b

where the extended system matrix is:

.. math::

   S =
   \begin{bmatrix}
   M & W_g^T & W_\gamma^T \\
   W_g & -\tfrac{C}{\theta \, dt^2} & 0 \\
   W_\gamma & 0 & -\tfrac{A}{\theta \, dt}
   \end{bmatrix}

and the right-hand side **b** encodes the current velocity, external forces,
carry-over multipliers, velocity source terms, and the Baumgarte position
correction.

When contacts are present, the contact impulses :math:`\Lambda_{NT}` enter through the velocity
correction:

.. math::

   u = S^{-1}(b + W^T \Lambda_{NT})

subject to the feasibility conditions on the contact impulses :math:`\Lambda_{NT}`, given by

.. math::
    \begin{aligned}
        &0 \le \gamma_{N}^+ \ \perp\ \Lambda_{N} \ge 0,\\
        &\lVert\boldsymbol{\Lambda}_{T,j}\rVert \le \mu_j\,\Lambda_{N,j},\\
        &\boldsymbol{\gamma}_{T,j}=\mathbf{0}\ \Rightarrow\ \lVert\boldsymbol{\Lambda}_{T,j}\rVert \le \mu_j\,\Lambda_{N,j}\qquad\text{(sticking)},\\
        &\boldsymbol{\gamma}^+_{T,j}\neq\mathbf{0}\ \Rightarrow\ \boldsymbol{\Lambda}_{T,j} = -\mu_j\,\Lambda_{N,j}\,\dfrac{\boldsymbol{\gamma}^+_{T,j}}{\lVert\boldsymbol{\gamma}^+_{T,j}\rVert}
        \quad\text{and}\quad \lVert\boldsymbol{\Lambda}_{T,j}\rVert = \mu_j\,\Lambda_{N,j}\qquad\text{(sliding)}.
    \end{aligned}

with post impact contact velocity :math:`\gamma_{NT}^{+}`.


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
   S is assembled once and factorised with a sparse LU. The free-motion velocity
   is the solution when no contacts are active. If there are no contacts the
   solver returns immediately.

**Jacobi sweep**
   Each iteration updates all impulse components simultaneously (in parallel,
   unlike Gauss-Seidel). The step size is scaled per row by the diagonal of the
   Delassus operator ``W M^-1 W^T``.

**Feasability projection**
   After each sweep, contact impulses are projected to satisfy:

   - **frictionless**: :math:`\Lambda_{N} \geq 0`
   - **frictional**: :math:`\|\Lambda_{T}\| \leq \mu\,\Lambda_{N}`

**Convergence**
   Checked every 10 iterations on the relative velocity residual. Controlled by
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
     - Convergence tolerances on the velocity residual
   * - ``solver.alpha``
     - Step-size multiplier for each Jacobi sweep (default 0.3)
   * - ``solver.nesterov``
     - Enable Nesterov momentum (default false)
   * - ``solver.warmstart``
     - Reuse previous-step impulses as initial guess (default true)
