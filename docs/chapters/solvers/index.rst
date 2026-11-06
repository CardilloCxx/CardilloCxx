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

.. toctree::
   :hidden:

   projected_jacobi
   projected_gauss_seidel
   condensed
   naive_schur
   interior_point
