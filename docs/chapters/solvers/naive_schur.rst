Naive Schur-Complement Solver (Reference)
=========================================

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

``solver.type = naive_schur`` (:cpp:class:`NaiveSchurSolver <cardillo::solver::NaiveSchurSolver>`)
eliminates every bilateral row -- springs, dampers and perfect constraints --
exactly by forming the dense Schur complement of the one-step problem onto
the contact rows, and then iterates on the contact percussions only. It is a
brute-force reference for the matrix-free solvers, in particular for
:doc:`condensed` with ``condensed.true_schur=true``, not a production solver.

Formulation
-----------

With the saddle-point matrix :math:`\mathcal{S}` of :doc:`projected_jacobi`
(factorized once per step by the same assembler) and the contact force
directions :math:`W_c` (columns), one solve with :math:`n_c` right-hand sides gives

.. math::

   X = \mathcal{S}^{-1} \begin{pmatrix} W_c \\ 0 \\ 0 \end{pmatrix},
   \qquad
   x = x^\mathrm{free} + X \pi_c ,

i.e. the new velocities and bilateral percussions are affine in the contact
percussions. The contact velocities are

.. math::

   \xi_c = c + G_{cc} \pi_c, \qquad
   G_{cc} = W_c^\top X_u = W_c^\top M^{-1} \big(W_c - W_b G_{bb}^{-1} G_{bc}\big) ,

with :math:`G_{bb} = W_b^\top M^{-1} W_b + \hat C_b`, :math:`G_{bc} = W_b^\top M^{-1} W_c` and :math:`c` the
contact velocity of the free solution (including restitution). The contact
law is solved by the projection
:math:`\pi_i = -\mathrm{proj}_{C_i}(\rho_i \xi_i - \pi_i)` per contact block,
with step lengths :math:`\rho_i` from the diagonal block :math:`G_{ii}` of the
Schur complement (``condensed.projection_rho_strategy``) times
``pj.alpha * pj.relaxation``, or by the local Alart--Curnier Newton solve on
:math:`G_{ii}` (``condensed.local_solve = newton``).

Difference to ``condensed.true_schur``
--------------------------------------

``condensed.true_schur`` never forms :math:`G_{cc}`. Each outer iteration first
solves the bilateral rows exactly for the current contact percussions, then
sweeps the contacts with the *mass-only* local blocks
:math:`W_i M^{-1} W_i^\top`. That is a block Gauss--Seidel iteration on the
full system with the bilateral rows as one super-block. Both iterations have
the same fixed points, but within a sweep a contact of the condensed solver
does not see the reaction of the bilateral rows to the other contacts, and
its step length ignores them. In the naive Schur solver, every local update
sees the exact bilateral response.

Configuration
-------------

No keys of its own; it reads the keys of the condensed and projected Jacobi
solvers so that both can be compared by switching only ``solver.type``. The contacts are swept
by block Gauss--Seidel:

- ``condensed.local_solve``, ``condensed.newton_*``,
  ``condensed.projection_rho_strategy``.
- ``pj.alpha``, ``pj.relaxation``, ``pj.tol_abs``, ``pj.tol_rel``,
  ``pj.max_iterations``, ``pj.warmstart``.
- ``pj.nesterov`` (and ``pj.nesterov_*``): Nesterov momentum with adaptive restart.

The convergence criterion is the same as for the condensed solver (scaled
change of the velocities between two sweeps).

Cost
----

One sparse LU of :math:`\mathcal{S}` (as projected Jacobi), :math:`n_c`
triangular solves, an :math:`n_v \times n_c` matrix :math:`X` and an
:math:`n_c \times n_c` matrix :math:`G_{cc}`; every sweep costs
:math:`\mathcal{O}(n_c^2)`. Use it for scenes with up to a few thousand
contact rows.
