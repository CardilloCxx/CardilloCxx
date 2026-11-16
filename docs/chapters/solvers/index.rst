Solvers
=======

Solving the contact problem is the research core of Cardillo. Each time step
produces a set of assembled operators from the collision and constraint pipeline
(see :doc:`../physics_pipeline`) and the solver's job is to find a new velocity
that satisfies the discrete Moreau time-stepping equations (see
:doc:`../moreau_time_stepping` for the derivation of the scaled system and RHS
terms such as the position-level force law and trajectory source velocities).
state that satisfies the equations of motion while respecting all contact and
constraint conditions simultaneously.

The contact problem is inherently non-smooth: contact forces only push (never
pull), and friction forces are bounded by a cone. Different solver architectures
handle this non-smooth feasibility problem in very different ways, which is why
the engine supports multiple backends.

All solvers share the same input from the assembler:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - **M**
     - Block-diagonal mass matrix (3x3 translational + 3x3 inertia per body)
   * - **u**
     - Stacked body velocities at the start of the step
   * - **f**
     - External and gyroscopic forces
   * - **Wg**, **Wgamma**
     - Constraint Jacobians (spring rows and damper rows)
   * - **C**, **A**
     - Compliance diagonals (spring and damper inverses)
   * - **W**
     - Contact Jacobian
   * - **mu**
     - Per-contact friction coefficients

And all solvers return the same output: the stacked body velocity vector **u**
for the new state, plus updated constraint multipliers **lambda_g** and
**lambda_gamma** for use in the next step.

Select a solver via ``solver.type`` in the config file (see :cpp:struct:`Config <cardillo::config::Config>`):

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Config value
     - Backend
     - Notes
   * - ``projected_jacobi`` (default)
     - :doc:`projected_jacobi` (:cpp:class:`ProjectedJacobiSolver <cardillo::solver::ProjectedJacobiSolver>`)
     - Production solver: direct factorization of the bilateral saddle-point
       system, projected iteration on the contacts. Most reliable.
   * - ``condensed``
     - :doc:`condensed` (:cpp:class:`CondensedSolver <cardillo::solver::CondensedSolver>`)
     - Production solver: block-sparse Schur complement; sweep strategy
       (Gauss-Seidel or graph-colored parallel), local solve (projection or
       semismooth Newton), Nesterov/Chebyshev acceleration and exact
       elimination of the bilateral rows (``condensed.true_schur``, needed
       for stiff rods).
   * - ``projected_gauss_seidel``
     - :doc:`projected_gauss_seidel` (:cpp:class:`ProjectedGaussSeidelSolver <cardillo::solver::ProjectedGaussSeidelSolver>`)
     - Small reference solver.
   * - ``naive_schur``
     - :doc:`naive_schur` (:cpp:class:`NaiveSchurSolver <cardillo::solver::NaiveSchurSolver>`)
     - Reference solver: dense Schur complement onto the contacts, all
       bilateral rows eliminated exactly. Small scenes only.
   * - ``qoco / clarabel / conicxx / scs``
     - :doc:`interior_point` (:cpp:class:`QocoSolver <cardillo::solver::QocoSolver>`, :cpp:class:`ClarabelSolver <cardillo::solver::ClarabelSolver>`, :cpp:class:`ConicxxSolver <cardillo::solver::ConicxxSolver>`, :cpp:class:`ScsSolver <cardillo::solver::ScsSolver>`)
     - Conic backends for the convexified contact law. Optional: each must be
       enabled at configure time with its CMake option
       (``CARDILLO_WITH_QOCO``, ``CARDILLO_WITH_CLARABEL``,
       ``CARDILLO_WITH_CONICXX``, ``CARDILLO_WITH_SCS``; all ``OFF`` by
       default). Selecting a backend that was not built throws a runtime
       error naming the option. ``conicxx`` additionally reuses its KKT
       factorization and warm-starts across steps when the active contact
       set is unchanged.
   * - ``mosek``
     - :doc:`interior_point` (``MosekSolver``)
     - Interior-point via the commercial MOSEK conic optimizer; same problem
       as ``clarabel``. Optional (``CARDILLO_WITH_MOSEK``), needs a MOSEK
       installation and license.

Accuracy measure: the natural map
---------------------------------

Accuracy is measured with the natural map of the contact law. For a contact
:math:`b` with percussion :math:`\lambda_b`, contact velocity :math:`\xi_b`
(restitution included) and a step :math:`\rho_b > 0`,

.. math::

   r_b = \rho_b^{-1}\Bigl(\lambda_b - \Pi_{K_\mu}\bigl(\lambda_b - \rho_b\,(\xi_b + s_b\,e_N)\bigr)\Bigr),

with the Coulomb cone :math:`K_\mu = \{\lVert\lambda_T\rVert \le \mu\lambda_N\}`
(frictionless contacts: :math:`\lambda_N \ge 0`) and the De Saxcé shift
:math:`s_b`:

- :math:`s_b = \mu\lVert\xi_{T}\rVert`: Coulomb's law, solved by all projection-based solvers.
- :math:`s_b = \mu\lVert\gamma_T^-\rVert`, frozen at the pre-step slip: the
  convexified law of the conic backends (:doc:`interior_point`).

:math:`r_b = 0` exactly when the contact satisfies its law, for any step and
also for non-unique (hyperstatic) percussions. :math:`r_b` is a velocity: for a
sticking contact it is the velocity error itself, whatever the step, for a
spurious percussion on a separating contact the velocity
:math:`\lambda_b/\rho_b` it would cause with the step :math:`\rho_b`. The
natural map is used in two places.

**Stopping test** (``pj.criterion = natural_map``, projected Jacobi). A
projection step of the solver, :math:`\lambda_b^+ = -\mathrm{proj}_C(\rho_b\xi_b - \lambda_b)`
with its own step :math:`\rho_b = \alpha R_b`, already evaluates the natural map
at its input in the Alart–Curnier form, which has the same zeros:
:math:`r_b = (\lambda_b - \lambda_b^+)/\rho_b` (both tangential rows need the
same step). The test costs a few operations per contact inside the projection
loop (``misc::naturalMapError2``) and is evaluated in every sweep. It scales
the normal and the tangential part separately, so that a normal error is not
hidden behind a large slip:

.. math::

   e = \max_b \max\!\left(
   \frac{|r_{b,N}|}{\varepsilon_\mathrm{abs} + \varepsilon_\mathrm{rel}|\xi_{b,N}|},
   \frac{\lVert r_{b,T}\rVert}{\varepsilon_\mathrm{abs} + \varepsilon_\mathrm{rel}\lVert\xi_{b,T}\rVert}
   \right) \le 1 ,

with ``pj.tol_abs`` and ``pj.tol_rel`` in velocity units.

**Final check** (``pj.final_check = natural_map``, every solver). After every
step, the time integrator evaluates the natural map of the stored percussions
at the returned velocities (``solver/natural_map.hpp``), with one step per
contact that does not depend on the solver, :math:`\rho_b = 1/\lambda_\mathrm{max}(G_b)`
from the contact's Delassus block :math:`G_b = W_b M^{-1} W_b^\top` (diagonal
mass); the conic backends have no step of their own. It reports raw residuals
in velocity units, independent of the tolerances: :math:`\max_b \max(|r_{b,N}|, \lVert r_{b,T}\rVert)`
for the solver's own law, the same for Coulomb's law, and the largest residual
:math:`W_i^\top u + b_i + c_i\Lambda_i` of the spring and damper rows. The
result is printed with ``debug.pj`` and available as
``PhysicsEngine::lastStepResidual()``. Near convergence the stopping test and
the final check agree; far from it the stopping test can be larger by up to
:math:`1/\alpha` (a spurious percussion counts with PJ's step
:math:`\alpha R_b`).

Comparing solvers:

- **Convergence** (work against accuracy): use the natural map of each
  solver's own law. It is zero at the solution of the solver's problem, so it
  measures how far each solver is from its solution with the same function.
  Coulomb's law is a harder problem than the convex convexified one.
- **Physics**: use the natural map of Coulomb's law. For the conic backends it
  contains the error of the convexified law, which does not decrease with the
  solver tolerance.

.. toctree::
   :hidden:

   projected_jacobi
   projected_gauss_seidel
   condensed
   naive_schur
   interior_point
