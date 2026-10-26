Interior-Point Solvers (QOCO, Clarabel and ConicXX)
====================================================

.. contents:: On this page
   :local:
   :depth: 2

QOCO, Clarabel and ConicXX formulate the one-step problem of the Moreau
scheme as a convex second-order cone program (SOCP) and solve it with a
primal-dual interior-point method. The bilateral force laws and frictionless
contacts are reproduced exactly. Coulomb friction is not a convex law, so the
conic backends solve a *convexified* friction law instead. Its error is
characterized below. Contact percussions always lie in the friction cone, but
the normal velocity of an active frictional contact is wrong by :math:`\mu` times
the change of its sliding speed within the step.

The derivation below follows the theory manual (``LaTeX/theory_doc.tex``,
section "Second-order cone formulation of the one-step problem"), which also
contains the proofs and the closed-form benchmark. Notation: time step
:math:`\Delta t`, percussions :math:`\pi = \Delta t\,\lambda` acting as
:math:`+W\pi`, scaled compliances
:math:`\hat C_g = \mathbf{C}/(\theta^2\Delta t^2)` (springs) and
:math:`\hat C_\gamma = \mathbf{A}/(\theta\Delta t)` (dampers), biases :math:`b_g`,
:math:`b_\gamma` of the position-level force law
(``DynamicsAssembler::springBias()``, ``damperBias()``), contact force directions
:math:`W_c`. Gyroscopic forces are explicit (:math:`f = f^\mathrm{ext} + f^\mathrm{gyr}(u_n)`),
see :ref:`ip-limitations`.

The one-step problem
--------------------

Unknowns are the new velocities :math:`u = u_{n+1}` and the percussions. The
discrete balance of momentum and the scaled force laws are

.. math::

   M u &= M u_n + \Delta t\,f + W_g \pi_g + W_\gamma \pi_\gamma + W_c \pi_c ,
   \\
   0 &= W_g^\top u + \hat C_g \pi_g + b_g ,
   \\
   0 &= W_\gamma^\top u + \hat C_\gamma \pi_\gamma + b_\gamma ,

and every closed contact satisfies the Signorini--Coulomb law between its
percussion :math:`\pi_i = (\pi_N, \pi_T)` and the velocity of Newton's impact
law

.. math::

   \xi_i = W_{c,i}^\top u + b_{c,i}, \qquad
   b_{c,i} = \chi_{c,i} + \mathrm{diag}(e_N, e_T, e_T)\,\gamma_i^-, \qquad
   \gamma_i^- = W_{c,i}^\top u_n + \chi_{c,i},

where :math:`\chi_c` is the velocity of static or driven partners and :math:`e_N`,
:math:`e_T` are the restitution coefficients.

Bilateral part: a convex QP
---------------------------

Without contacts, the one-step problem is the optimality system of the convex
quadratic program

.. math::

   \min_{u,\,\pi_g,\,\pi_\gamma}\;
   \tfrac{1}{2} u^\top M u - u^\top (M u_n + \Delta t\,f)
   + \tfrac{1}{2} \pi_g^\top \hat C_g \pi_g + \tfrac{1}{2} \pi_\gamma^\top \hat C_\gamma \pi_\gamma
   \quad\text{s.t.}\quad
   W_g^\top u + \hat C_g \pi_g + b_g = 0, \quad
   W_\gamma^\top u + \hat C_\gamma \pi_\gamma + b_\gamma = 0 .

The objective is Gauss' principle for the velocity jump,
:math:`\tfrac12\|u - u^\mathrm{free}\|_M^2`, plus the scaled complementary energy
of the compliant elements. With multipliers :math:`y_g`, :math:`y_\gamma` of the
equality constraints, the stationarity conditions are

.. math::

   M u - (M u_n + \Delta t\,f) + W_g y_g + W_\gamma y_\gamma &= 0 ,
   \\
   \hat C_g \pi_g + \hat C_g^\top y_g &= 0 ,
   \\
   \hat C_\gamma \pi_\gamma + \hat C_\gamma^\top y_\gamma &= 0 .

The transposes come from the transposed constraint matrix. Since the scaled
compliances are symmetric (diagonal in the code), the last two rows read
:math:`\hat C(\pi + y) = 0`. Together with the constraints, :math:`u` and the
percussions :math:`-y` therefore solve the momentum balance and the force laws
above. For compliant rows, :math:`\pi = -y`. For perfect constraints
(:math:`\hat C = 0`), the primal :math:`\pi` does not enter the problem and is not
determined; the percussion is :math:`-y` only.

Keeping :math:`\pi` as a primal unknown, instead of eliminating it into a
penalty term :math:`\tfrac12\|W^\top u + b\|^2_{\hat C^{-1}}`, keeps the problem
well defined in the rigid limit :math:`\hat C \to 0`. The objective needs a
symmetric positive definite mass matrix, which is why implicit gyroscopic
forces (non-symmetric :math:`M_\mathrm{eff}`) are not supported.

Coulomb friction as a cone complementarity problem
--------------------------------------------------

For a contact with :math:`\mu > 0`, let
:math:`K_\mu = \{\pi : \|\pi_T\| \le \mu\,\pi_N\}` be the Coulomb cone and
:math:`K_\mu^* = \{\xi : \mu\|\xi_T\| \le \xi_N\}` its dual cone. De Saxcé's
reformulation (De Saxcé & Feng 1998) states that the Signorini--Coulomb law is
equivalent to

.. math::

   K_\mu^* \ni \hat\xi \;\perp\; \pi \in K_\mu,
   \qquad
   \hat\xi = \xi + \begin{pmatrix} \mu \|\xi_T\| \\ 0 \\ 0 \end{pmatrix} .

The cases are separation (:math:`\pi = 0`, :math:`\xi_N \ge 0`), sticking
(:math:`\pi` interior, :math:`\hat\xi = 0`) and sliding (:math:`\pi` on the cone
boundary, :math:`\xi_N = 0`, :math:`\xi_T` opposite to :math:`\pi_T`). The shift
:math:`\mu\|\xi_T\|` is not affine in :math:`u`, so the problem is not the
optimality system of a convex program. Without the shift, the conditions
describe Anitescu's associated relaxation, which requires the separation
velocity :math:`\xi_N = \mu\|\xi_T\|` while sliding.

Convexification
---------------

The conic backends freeze the shift at the pre-step velocity,
:math:`c_i = \mu_i \|\gamma_{T,i}^-\|`. The shifted velocity
:math:`\xi_i + (c_i, 0, 0)^\top` is then affine in :math:`u`, and the one-step
problem becomes the SOCP

.. math::

   \min_{u,\,\pi_g,\,\pi_\gamma}\; \text{(objective above)}
   \quad\text{s.t.}\quad
   \text{(bilateral equalities)}, \quad
   W_{c,i}^\top u + b_{c,i} + (c_i, 0, 0)^\top \in K_{\mu_i}^* , \quad
   w_{N,j}^\top u + b_{N,j} \ge 0 ,

with one cone constraint per frictional contact :math:`i` and one inequality per
frictionless contact :math:`j`. The multipliers of the cone constraints are
exactly the contact percussions :math:`\pi_{c,i} \in K_{\mu_i}`, and the KKT
conditions are the balance of momentum, the bilateral force laws, Signorini's
law for frictionless contacts and the *convexified* friction law

.. math::

   K_\mu^* \ni \xi + (c, 0, 0)^\top \;\perp\; \pi \in K_\mu .

Scaling to the standard cone
----------------------------

The backends work with the self-dual standard cone
:math:`L = \{(a, b) \in \mathbb{R}\times\mathbb{R}^2 : \|b\| \le a\}`. With the
per-contact scaling :math:`S_\mu = \mathrm{diag}(1, \mu, \mu)`,

.. math::

   \pi \in K_\mu \iff z = S_\mu^{-1}\pi \in L,
   \qquad
   \hat\xi \in K_\mu^* \iff s = S_\mu \hat\xi \in L,
   \qquad
   s^\top z = \hat\xi^\top \pi ,

so that the convexified law becomes the standard cone complementarity
:math:`L \ni s \perp z \in L`. The normal row keeps its physical scale and the
scaling stays bounded for :math:`\mu \to 0`, unlike the equivalent
:math:`\mathrm{diag}(1/\mu, 1, 1)`, which blows up the normal row. Since the
normal entry of :math:`S_\mu` is one, the shift is not scaled. Frictionless
contacts are single rows in the non-negative orthant.

Standard form
-------------

With :math:`x = (u, \pi_g, \pi_\gamma)`, :math:`S = \mathrm{blkdiag}(S_{\mu_i})`
(one for frictionless rows) and :math:`\bar c = ((c_i, 0, 0))_i`, the assemblers
(``*_assembler.cpp``, scaling and shift in ``conic_scaling.hpp``) build the
matrices below. The sans-serif names are those of the backend interfaces, unrelated
to the damper compliance :math:`\mathbf{A}` and the time step.

.. math::

   \mathsf{P} = \begin{pmatrix}
     M & 0 & 0 \\
     0 & \hat C_g & 0 \\
     0 & 0 & \hat C_\gamma
   \end{pmatrix},
   \qquad
   \mathsf{q} = -\begin{pmatrix}
     M u_n + \Delta t\,f \\ 0 \\ 0
   \end{pmatrix},

.. math::

   \mathsf{A}_\mathrm{eq} = \begin{pmatrix}
     W_g^\top & \hat C_g & 0 \\
     W_\gamma^\top & 0 & \hat C_\gamma
   \end{pmatrix},
   \qquad
   \mathsf{b}_\mathrm{eq} = -\begin{pmatrix} b_g \\ b_\gamma \end{pmatrix},
   \qquad
   \mathsf{G} = \begin{pmatrix} -S W_c^\top & 0 & 0 \end{pmatrix},
   \qquad
   \mathsf{h} = S b_c + \bar c .

The cone slack is :math:`s = \mathsf{h} - \mathsf{G}x = S(W_c^\top u + b_c) + \bar c = S\hat\xi`,
and the cone is :math:`\mathcal{K} = \mathbb{R}_{\ge 0}^{n_N} \times L \times \dots \times L`.
QOCO solves

.. math::

   \min_x\; \tfrac{1}{2} x^\top \mathsf{P} x + \mathsf{q}^\top x
   \quad\text{s.t.}\quad
   \mathsf{A}_\mathrm{eq} x = \mathsf{b}_\mathrm{eq}, \quad
   \mathsf{h} - \mathsf{G} x = s \in \mathcal{K} .

Clarabel and ConicXX stack the equality rows as a zero cone on top:
:math:`\mathsf{A} = (\mathsf{A}_\mathrm{eq}; \mathsf{G})`,
:math:`\mathsf{b} = (\mathsf{b}_\mathrm{eq}; \mathsf{h})`, cone
:math:`\{0\}^{m_g + m_\gamma} \times \mathcal{K}`.

KKT conditions and recovery of the solution
-------------------------------------------

With multipliers :math:`y = (y_g, y_\gamma)` of the equality rows and :math:`z`
of the cone rows (for Clarabel and ConicXX, :math:`y` is the leading zero-cone
part of their dual vector), the optimal primal-dual solution satisfies

.. math::

   \mathsf{P}x + \mathsf{q} + \mathsf{A}_\mathrm{eq}^\top y + \mathsf{G}^\top z &= 0
   \quad\text{(stationarity)}
   \\
   \mathsf{A}_\mathrm{eq} x &= \mathsf{b}_\mathrm{eq}
   \quad\text{(force laws)}
   \\
   \mathsf{h} - \mathsf{G} x &= s
   \quad\text{(cone slack)}
   \\
   \mathcal{K} \ni s &\perp z \in \mathcal{K}
   \quad\text{(complementarity)}

Expanding the stationarity condition with the concrete matrices gives

.. math::

   M u - M u_n - \Delta t\,f + W_g\,y_g + W_\gamma\,y_\gamma - W_c S z &= 0 ,
   \\
   \hat C_g\,\pi_g + \hat C_g^\top y_g &= 0 ,
   \\
   \hat C_\gamma\,\pi_\gamma + \hat C_\gamma^\top y_\gamma &= 0 .

The solution of the one-step problem with the convexified friction law is
therefore

.. math::

   u_{n+1} = u, \qquad
   \pi_c = S z, \qquad
   \pi_g = -y_g, \qquad
   \pi_\gamma = -y_\gamma ,

and the first row is the discrete balance of momentum
:math:`M u_{n+1} = M u_n + \Delta t\,f + W_g\pi_g + W_\gamma\pi_\gamma + W_c\pi_c`.
The solvers return :math:`u`, store :math:`S z` as ``Contact::last_impulse`` and
:math:`-y` as bilateral percussions (``Lambda_g``, ``Lambda_gamma``), which also
covers perfect constraints.

Error of the convexified friction law
-------------------------------------

Every active frictional contact (:math:`\pi \ne 0`) of the conic solvers
satisfies

.. math::

   \xi_N = \mu \bigl(\|\xi_T\| - \|\gamma_T^-\|\bigr)

instead of Signorini's :math:`\xi_N = 0`. A separated contact satisfies
:math:`\xi_N \ge \mu(\|\xi_T\| - \|\gamma_T^-\|)` instead of :math:`\xi_N \ge 0`.
The friction cone, the direction of the friction percussion and the
stick/slip alternative are those of Coulomb's law. Consequences:

- **Exact** for permanent sticking and steady sliding (the sliding speed does
  not change within the step).
- **Persistent contact:** a decelerating slip pushes the bodies into each
  other and an accelerating slip lifts them off. The normal velocity error is
  :math:`\mathcal{O}(\Delta t)` per step, and the gap error telescopes to
  :math:`\mu\,\Delta t\,(\|\gamma_{T,\mathrm{end}}\| - \|\gamma_{T,\mathrm{start}}\|)`,
  i.e., :math:`\mathcal{O}(\Delta t)` without accumulation.
- **Frictional impacts:** the sliding speed jumps within one step, and the
  post-impact normal velocity has an :math:`\mathcal{O}(1)` error that does not
  vanish for :math:`\Delta t \to 0`.

The theory manual quantifies this with a solid sphere on a plane
(:math:`\kappa = 1 + mR^2/\Theta = 7/2`), checked against all three
backends in ``tests/test_convexified_friction.cpp``:

- A ball **sliding** with :math:`v_0` and no spin sinks with
  :math:`\kappa\mu^2 g\,\Delta t` and ends up :math:`\mu\,\Delta t\,v_0` deep in the
  plane (0.6 mm for :math:`\mu = 0.3`, :math:`v_0 = 2\,\mathrm{m/s}`,
  :math:`\Delta t = 1\,\mathrm{ms}`). Its rolling velocity is correct.
- A **spinning ball dropped** onto the plane rebounds with the effective
  restitution :math:`e^* = (e - \kappa\mu^2)/(1 + \kappa\mu^2)` (sliding impact)
  or :math:`e^* = e - \mu\,\gamma_T^-/V` (sticking impact) instead of :math:`e`,
  independently of :math:`\Delta t`. For :math:`\mu = 0.3` and :math:`e = 0.5`,
  :math:`e^* = 0.14`: the ball reaches 8 % of the correct rebound height. For
  :math:`e < \kappa\mu^2`, :math:`e^* < 0`, and the ball moves into the plane
  after the impact. After an inelastic impact (:math:`e = 0`), it sinks by
  :math:`\mu\,\Delta t\,\gamma_T^-` until it rolls (1.5 mm for
  :math:`\gamma_T^- = 5\,\mathrm{m/s}`, :math:`\Delta t = 1\,\mathrm{ms}`): the
  velocity error of the impact step is :math:`\mathcal{O}(1)`, but the
  penetration is :math:`\mathcal{O}(\Delta t)`.

Use projected Jacobi or the condensed solver, which solve the exact law, when
frictional impacts matter.

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

.. _ip-limitations:

Limitations
-----------

- **Convexified friction.** Frictional contacts follow the convexified law
  above, not Coulomb's law; frictional impacts have an error that does not
  vanish with the step size.
- **QOCO/Clarabel: no warmstart, rebuilt every step.** The last result might
  not be feasible in the next timestep, so it cannot directly be used as
  warmstart, and when contacts appear or disappear the sparsity pattern of
  :math:`\mathsf{G}`/:math:`\mathsf{A}` changes, requiring a full solver re-setup every step.
  ConicXX does not have this limitation in the common case where the contact
  set is momentarily stable (see above).
- **No implicit gyroscopy** (``moreau.implicit_gyroscopy``). This setting is
  silently ignored by all three interior-point backends (the QP needs a
  symmetric mass matrix); use PJ or the condensed solver if you need it.

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
