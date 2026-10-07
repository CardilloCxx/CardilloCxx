SCS (ADMM Conic Solver)
=======================

.. contents:: On this page
   :local:

``solver.type = scs`` (:cpp:class:`ScsSolver <cardillo::solver::ScsSolver>`)
solves exactly the cone program of the interior-point backends (see
:doc:`interior_point`) with `SCS <https://github.com/cvxgrp/scs>`_, a
first-order solver that applies Douglas-Rachford splitting (ADMM) to the
homogeneous self-dual embedding, with Anderson acceleration and adaptive step
scaling. It is only compiled with ``-DCARDILLO_WITH_SCS=ON``; SCS (pinned to
3.3.1) is fetched and built from source.

Problem data
------------

The data comes unchanged from the Clarabel assembler: :math:`x = [u;
\lambda_g; \lambda_\gamma]`, a diagonal :math:`P`, the P\ :sub:`μ`-scaled
contact rows and the friction shift in :math:`b`. SCS's standard form

.. math::

   \min_x \tfrac12 x^\top P x + q^\top x
   \quad \text{s.t.} \quad Ax + s = b,\; s \in \{0\}^p \times \mathbb{R}_+^l
   \times L^3 \times \dots \times L^3

is Clarabel's, and SCS's mandatory cone order (zero, nonnegative, second-order,
scalar component first) is the assembler's row order, so nothing is permuted.
The dual :math:`y` satisfies :math:`Px + q + A^\top y = 0`; contact impulses are
:math:`y \circ S_\mu`, as for the other conic backends.

Only SCS's direct linear-system backend is used: the quasi-definite KKT matrix
is ordered with AMD and factorized with QDLDL in ``scs_init``. Because
:math:`A` changes with the configuration, every step runs a fresh
``scs_init`` (timer ``SCS Setup``). Every adaptive scale update inside the
iteration triggers a numerical refactorization of the same cost, which SCS does
not include in its own ``lin_sys_time``.

Warm start
----------

Unlike the interior-point backends, SCS can be warm started (``scs.warm_start``,
default ``on``):

* :math:`x_0 = [u_{\text{old}}; \lambda_g; \lambda_\gamma]` (start-of-step
  velocities and the stored compliant multipliers);
* :math:`y_0`: stationarity in the compliant columns gives :math:`y = \lambda`
  for the spring/damper rows; contact rows use the impulses tracked across steps
  by the contact matcher (the same storage PGS/PJ warm start from), divided by
  :math:`S_\mu` and projected onto the cones;
* :math:`s_0 = \Pi_K(b - A x_0)`, zero on the equality rows.

With ``scs.carry_scale`` (default ``on``) each step starts from the previous
step's final adapted dual scale.

Build notes
-----------

SCS vendors AMD and QDLDL under their original symbol names, which already exist
in the process (ConicXX's QDLDL in ``libcardillo``, QOCO's AMD in ``libqoco``).
``libscsdir`` is therefore a shared library whose exports are restricted to the
public ``scs_*`` API (``cmake/scs_exports.map``), so each library binds to its
own copy. ``scs_int``/``scs_float`` are checked against Eigen's index type and
``real_t`` at compile time.

SCS compiles Anderson acceleration only with BLAS/LAPACK, which is therefore
enabled when found (a warning is printed otherwise). SCS is built without
OpenMP unless ``-DCARDILLO_SCS_OPENMP=ON``: its OpenMP loops over the
three-dimensional contact cones fork a full thread team every iteration and made
small scenes several times slower.

Config keys
-----------

Negative numeric values (the default) keep SCS's own defaults.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Key
     - Meaning
   * - ``scs.eps_abs``, ``scs.eps_rel``
     - Stopping tolerances; default to ``pj.tol_abs`` / ``pj.tol_rel``. Not
       comparable one-to-one with interior-point tolerances.
   * - ``scs.eps_infeas``, ``scs.max_iters``, ``scs.time_limit_secs``
     - SCS settings of the same name
   * - ``scs.alpha``, ``scs.scale``, ``scs.rho_x``
     - Relaxation, initial dual scale, primal scale
   * - ``scs.adaptive_scale``, ``scs.normalize``
     - Adaptive scaling and data equilibration (both default ``on``)
   * - ``scs.acceleration_lookback``, ``scs.acceleration_interval``
     - Anderson acceleration memory (``0`` disables) and interval
   * - ``scs.warm_start``, ``scs.carry_scale``
     - See above (both default ``on``)
   * - ``scs.stats_csv``
     - Per-step statistics (iterations, setup/solve/linear-system/cone/AA
       times, accepted and rejected AA steps, scale updates, residuals)
   * - ``scs.dump_dir``, ``scs.dump_every``
     - Write every N-th step's cone program and warm start for offline replay
       with ``build/tests/conic_replay`` (see ``tools/conic_backend_benchmark.py``)

A status of *solved inaccurate* is accepted and counted (summary printed at
exit); any other non-solved status throws, as for the other conic backends.
``debug.pj = true`` enables SCS's own log.
