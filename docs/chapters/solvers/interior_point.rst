Interior-Point Solvers (QOCO, Clarabel and ConicXX)
====================================================

.. contents:: On this page
   :local:
   :depth: 2

QOCO, Clarabel and ConicXX all formulate the contact and constraint problem as
a convex second-order cone program (SOCP) and solve it with a primal-dual
interior-point method. Unlike the iterative projection solvers (PJ, PGS), they
enforce contact feasibility as a mathematical constraint rather than a
post-projection step, which means they can satisfy the Coulomb friction cone
exactly.

The QP formulation
------------------

The Moreau step (velocity :math:`x = (u,\,\Lambda_g,\,\Lambda_\gamma)`,
scaled multipliers as defined in :doc:`projected_gauss_seidel`) is cast as
a convex quadratic program:

.. math::

   \min_{x}\; \tfrac{1}{2}\,x^\top P\,x + c^\top x

subject to equality constraints from the compliant force laws and cone
constraints from contact.

**Objective** (:math:`P`, :math:`c`)

.. math::

   P = \begin{pmatrix}
     M_\mathrm{eff} & 0 & 0 \\
     0 & \tfrac{C}{\theta\,dt^2} & 0 \\
     0 & 0 & \tfrac{A}{\theta\,dt}
   \end{pmatrix},
   \qquad
   c = -\begin{pmatrix}
     M u_n + dt\,f^\mathrm{ext} \\[2pt] 0 \\ 0
   \end{pmatrix}

**Constraint equations** (:math:`A`, :math:`b`) -- spring and damper rows

.. math::

   A = \begin{pmatrix}
     W_g^\top & -\tfrac{C}{\theta\,dt^2} & 0 \\[2pt]
     W_\gamma^\top & 0 & -\tfrac{A}{\theta\,dt}
   \end{pmatrix},
   \qquad
   b = \begin{pmatrix}
     -\tfrac{C\,\Lambda_{g,n}}{\theta\,dt^2}
     - \tfrac{1-\theta}{\theta}W_g^\top u_n
     - \tfrac{1}{\theta}v_g \\[4pt]
     - \tfrac{1-\theta}{\theta}W_\gamma^\top u_n
     - \tfrac{1}{\theta}v_\gamma
   \end{pmatrix}

These enforce the compliant constraint laws :math:`Ax = b` as hard equalities
(hard constraints simply have :math:`C = 0` or :math:`A = 0`, making those rows
purely kinematic).

The Coulomb friction cone
-------------------------

For a single frictional contact :math:`i` with friction coefficient
:math:`\mu_i`, the admissible impulse set is the scaled second-order cone:

.. math::

   K_{\mu_i} = \bigl\{\,\Lambda_{NT} \in \mathbb{R}^3 \;\big|\;
   \|\Lambda_{T}\| \leq \mu_i\,\Lambda_{N}\,\bigr\}

Enforcing :math:`\Lambda_{NT} \in K_{\mu_i}` directly requires a scaled cone, which
complicates the interior-point system. Following De Saxce and Acary, the
implementation introduces the per-contact scaling matrix:

.. math::

   S_{\mu_i} = \mathrm{diag}\!\left(\tfrac{1}{\mu_i},\, 1,\, 1\right)

and the change of variables :math:`z_i = S_{\mu_i}^{-1}\Lambda_{NT, i}`, so that the
scaled impulse :math:`z_i` lives in the standard SOC :math:`L`:

.. math::

   \Lambda_{NT,i} \in K_{\mu_i}
   \quad\Longleftrightarrow\quad
   z_i = S_{\mu_i}^{-1}\Lambda_{NT, i} \in L

On the velocity side the contact velocity :math:`\gamma_i = W_{NT,i}^\top u`
is corrected by the De Saxce term :math:`\bar{\gamma}_i` which is approximated from velocity of last timestep as
:math:`\bar{\gamma} = \left( \mu\|\gamma_T\|, 0, 0 \right)^\top`.

This allows the frictional contact condition to be written as the standard cone
complementarity :math:`L \ni s_i \perp z_i \in L`.
Frictionless contacts reduce to a single row in the non-negative orthant.

**Contact friction**

The contact friction conditions are expressed as the standard cone complementarity
:math:`h_i - G_i\,x =s_i \in L` for each frictional contact :math:`i`.

.. math::

   G = \begin{pmatrix} -S_{\mu}^\top W_{NT}^\top &  0 & 0  \end{pmatrix}
   \qquad
   h = S_{\mu}^\top\!\left(\gamma^\mathrm{ext} + \bar{\gamma}\right)


where :math:`\gamma^\mathrm{ext}` is the contact velocity contribution
from non-physical moving bodies.

Because the non-convex coulomb friction problem is solved as a convexified problem,
the interior-point method does not solve the exact coloumb friction problem.

KKT conditions
--------------

The optimal primal-dual solution
:math:`(x^*,\, y^*,\, s^*,\, z^*)` satisfies:

.. math::

   P\,x + c + A^\top y + G^\top z &= 0
   \quad\text{(stationarity)}
   \\
   A\,x &= b
   \quad\text{(primal: constraint equations)}
   \\
   h - G\,x &= s
   \quad\text{(primal: cone slack)}
   \\
   L \ni s^i &\perp z^i \in L \quad \forall\, i
   \quad\text{(complementarity)}

Expanding the stationarity condition with the concrete matrices gives:

.. math::

   M_\mathrm{eff}\,u + W_g\,y_g + W_\gamma\,y_\gamma - W_{NT} S_\mu\, z
   &= M u_n + dt\,f^\mathrm{ext}
   \\
   \tfrac{C}{\theta\,dt^2}\,\Lambda_g - \tfrac{C}{\theta\,dt^2}\,y_g &= 0
   \;\Rightarrow\; y_g = \Lambda_g
   \\
   \tfrac{A}{\theta\,dt}\,\Lambda_\gamma - \tfrac{A}{\theta\,dt}\,y_\gamma &= 0
   \;\Rightarrow\; y_\gamma = \Lambda_\gamma

So the dual variables :math:`y_g` and :math:`y_\gamma` coincide with the
primal multipliers, and the first KKT equation recovers the discrete momentum
balance. The constraint rows enforce the compliant force laws, and the
complementarity condition enforces Signorini + Coulomb.

QOCO vs Clarabel vs ConicXX
---------------------------

All three backends solve the same physical problem with the same matrices
(ConicXX and Clarabel share the exact same standard form). The differences
are in standard form, implementation, and build requirements.

.. list-table::
   :header-rows: 1
   :widths: 20 30 30 30

   * -
     - QOCO
     - Clarabel
     - ConicXX
   * - Standard form
     - :math:`\min \tfrac{1}{2}x^\top Px + c^\top x` s.t. :math:`Ax=b`, :math:`Gx \in \mathcal{K}`
     - :math:`\min \tfrac{1}{2}x^\top Px + q^\top x` s.t. :math:`Ax \in \mathcal{K}`
     - :math:`\min \tfrac{1}{2}x^\top Px + q^\top x` s.t. :math:`Ax + s = b,\ s \in \mathcal{K}` (same layout as Clarabel)
   * - Equality handling
     - Separate :math:`A`, :math:`b` matrices
     - Zero-cone rows prepended to :math:`A`
     - Zero-cone rows prepended to :math:`A` (``ConeSpec::zero_dim``)
   * - Iterative refinement
     - Not supported
     - ``solver.iter_ref_iters`` refinement passes per KKT solve
     - ``solver.iter_ref_iters`` refinement passes per KKT solve
   * - Algorithm
     - Primal-dual interior point (homogeneous embedding)
     - Homogeneous self-dual interior point
     - Homogeneous self-dual interior point (Nesterov-Todd scaling, Mehrotra
       predictor-corrector)
   * - Warmstart
     - None (restarts from central path each step)
     - None
     - **Yes** -- ``Solver::updateData()`` reuses the KKT sparsity pattern
       and factorization symbolic analysis, and ``Settings::warm_start``
       seeds the next solve from the previous iterate, whenever the active
       contact set hasn't changed since the last step

QOCO backend selection
----------------------

QOCO resolves its backend library at runtime:

- ``qoco.backend = auto`` -- tries CUDA first, falls back to CPU
- ``qoco.backend = cpu`` -- CPU only
- ``qoco.backend = cuda`` -- CUDA only

ConicXX warm start and incremental reuse
-----------------------------------------

Unlike QOCO and Clarabel, :cpp:class:`ConicxxSolver <cardillo::solver::ConicxxSolver>`
keeps a single ``conicxx::Solver`` instance alive for the whole simulation.
Every step it rebuilds ``P``/``q``/``A``/``b``/cone spec from the assembler as
usual, but instead of unconditionally re-creating the solver, it first tries
``conicxx::Solver::updateData()`` (a cheap values-only update that reuses the
existing KKT sparsity pattern and factorization symbolic analysis). Only when
that fails -- i.e. the active contact set changed and the sparsity pattern of
``A`` (or the cone composition; see below) no longer matches -- does it fall
back to a full ``conicxx::Solver::setup()``, exactly mirroring the usage
pattern documented in ConicXX's own ``solver.h``. When ``conicxx.warm_start``
is enabled, ConicXX also seeds each solve from its own previous iterate
automatically (no explicit bookkeeping needed on the CardilloMPI side).

Note that ``updateData()`` only checks the sparsity pattern of ``A``, not the
cone composition itself -- a contact whose friction coefficient crosses zero
between steps can switch a row between the nonnegative orthant and a
second-order cone block while coincidentally leaving ``A``'s row count and
nonzero pattern unchanged. ``ConicxxSolver`` guards against this explicitly by
also comparing the freshly computed ``ConeSpec`` against the one used at the
last ``setup()``, and forcing a full re-setup if it differs.

Limitations
-----------

- **QOCO/Clarabel: no warmstart, rebuilt every step.** The last result might
  not be feasible in the next timestep, so it cannot directly be used as
  warmstart, and when contacts appear or disappear the sparsity pattern of
  :math:`G`/:math:`A` changes, requiring a full solver re-setup every step.
  ConicXX does not have this limitation in the common case where the contact
  set is momentarily stable (see above).
- **No implicit gyroscopy** (``moreau.implicit_gyroscopy``) and no
  lambda-theta integration (``moreau.lambda_theta``). These settings are
  silently ignored by all three interior-point backends; use PJ or PGS if you
  need them.

Config keys
-----------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Key
     - Effect
   * - ``solver.type = qoco`` / ``clarabel`` / ``conicxx``
     - Select backend
   * - ``qoco.backend``
     - ``auto`` / ``cpu`` / ``cuda`` (QOCO only)
   * - ``solver.max_iterations``
     - Maximum interior-point iterations
   * - ``solver.tol_abs``, ``solver.tol_rel``
     - Primal/dual feasibility and gap tolerances
   * - ``solver.kkt_static_reg``
     - Static KKT regularisation constant (all three)
   * - ``solver.kkt_dynamic_reg``
     - Dynamic KKT regularisation delta (all three)
   * - ``solver.iter_ref_iters``
     - Iterative refinement passes per KKT solve (Clarabel, ConicXX)
   * - ``conicxx.warm_start``
     - Enable/disable ConicXX's cross-step warm start and KKT reuse
       (ConicXX only, default ``on``)
   * - ``conicxx.linear_solver``
     - Which sparse LDL^T backend factorizes the KKT system: ``eigen``
       (``Eigen::SimplicialLDLT``, always available), ``qdldl`` (the backend
       QOCO/Clarabel use -- ConicXX's own default), or ``regularized_ldlt``
       (a custom Davis/ECOS-style backend correcting bad pivots inline
       during elimination). ConicXX only; ``qdldl``/``regularized_ldlt``
       silently fall back to ``eigen`` if ConicXX wasn't built with
       ``CONICXX_WITH_QDLDL``.

       .. warning::
          ``regularized_ldlt`` was found to need dramatically more
          interior-point iterations than ``eigen``/``qdldl`` on large,
          heavily-contacted scenes (a domino-toppling scene with ~13,600
          simultaneous frictional contacts didn't converge within 1200+
          iterations, where ``eigen``/``qdldl`` converge in ~80) despite
          identical per-factorization cost -- likely a convergence-quality
          issue in how it handles many simultaneous near-degenerate pivots,
          not a raw performance regression. Prefer ``qdldl`` (the default)
          or ``eigen`` until this is resolved upstream.
   * - ``conicxx.dynamic_reg_eps``
     - Pivot-magnitude threshold triggering a dynamic regularization bump
       (ConicXX only)
   * - ``conicxx.refine_tol``
     - Target relative residual for iterative refinement (ConicXX only)
   * - ``conicxx.max_step_fraction``
     - Fraction-to-boundary step-length safety factor (ConicXX only)
   * - ``conicxx.equilibrate``
     - Enable/disable Ruiz equilibration (ConicXX only, default ``on``)
   * - ``conicxx.equilibrate_max_iter``, ``conicxx.equilibrate_min_scale``,
       ``conicxx.equilibrate_max_scale``
     - Ruiz equilibration tuning knobs (ConicXX only)
   * - ``conicxx.validate_inputs``
     - Enable/disable ConicXX's input validation (ConicXX only, default
       ``on``)
