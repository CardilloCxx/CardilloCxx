Projected Gauss-Seidel
======================

.. contents:: On this page
   :local:
   :depth: 2

The core problem
----------------

The Moreau time-stepping scheme produces the scaled linear system used by
the solvers. In compact form:

.. math::

  \mathcal{S}\,\mathbf{x}_{n+1} = \mathbf{b}_n

See the Moreau chapter :doc:`../moreau_time_stepping` for the full derivation and the detailed block
definitions of :math:`\mathcal{S}` and :math:`\mathbf{b}_n`.

Condensed formulation
---------------------

PGS iterates on the condensed (impulse) form of the one-step problem, in the
notation of the theory manual (``LaTeX/theory_doc.tex``, Section "The condensed
form of the one-step problem"). All block percussions act on the velocities with
the force directions of ``DynamicsAssembler`` as they are,

.. math::

   \pi = \begin{pmatrix} \pi_g \\ \pi_\gamma \\ \pi_c \end{pmatrix},
   \qquad
   W = \begin{pmatrix} W_g & W_\gamma & W_c \end{pmatrix},
   \qquad
   u_{n+1} = u_\mathrm{free} + \delta u, \qquad \delta u = M^{-1} W \pi,

with the bilateral percussions :math:`\pi_g = h\lambda_g`,
:math:`\pi_\gamma = h\lambda_\gamma` and the contact percussions
:math:`\pi_c = (\pi_N, \pi_T)`, :math:`\pi_N \geq 0` (``PgsAssembler::W()``).
The row velocities are

.. math::

   \xi = \xi^\mathrm{free} + W^\top \delta u + \hat{C}\pi,

with :math:`\xi^\mathrm{free} = W^\top u_\mathrm{free} + b` (``PgsAssembler::rhs()``)
and :math:`\hat C = \mathrm{diag}(C/(\theta^2 h^2),\, A/(\theta h),\, 0)`
(``PgsAssembler::C()``). For a bilateral block, :math:`\xi_i = 0` is the
discrete force law; for a contact block, :math:`\xi_i` is the contact
velocity. This avoids ever forming the Delassus operator explicitly.

The block diagonal preconditioner
----------------------------------

PGS preconditions each block :math:`j` with the inverse of the diagonal of
:math:`S_{jj} = W_j^\top M^{-1} W_j + \hat{C}_j` (``PgsAssembler::DinvDiag()``).

The PGS algorithm
-----------------

PGS differs from Projected Jacobi in that it updates each block sequentially
and propagates the velocity correction :math:`\delta u` immediately,
so later blocks in the same sweep already see the effect of earlier updates.
Contact blocks use the negative variant of the projection,
:math:`\pi_N = -\mathrm{proj}_{\mathbb{R}_{\leq 0}}(\rho_N\xi_N - \pi_N)` and
:math:`\pi_T = -\mathrm{proj}_{D(\mu\pi_N)}(\rho_T\xi_T - \pi_T)`
(``misc/contact_projection.hpp``, shared by all projection-based solvers).

.. code-block:: text

   Algorithm: Projected Gauss-Seidel
   ---------------------------------
   Input : rhs = xi_free, W, M^{-1}, D^{-1}, k_max, relaxation
   Output: u_{n+1}

   1. Initialize pi <- warm start (pi_g, pi_gamma, contact percussions)
   2. du <- M^{-1} W pi

   3. for k = 0 to k_max do
   4.     for each block j in order do
   5.         xi_j   <- rhs_j + W_j^T du + C_j pi_j
   6.         step_j <- omega * D_j^{-1} * xi_j
   7.         if j is contact block then step_j <- alpha * step_j

   8.         if j is bilateral then
   9.             pi_j_new <- pi_j - step_j
   10.        else
   11.            pi_j_new <- -proj_{C_j}(step_j - pi_j)

   12.        du  <- du + M^{-1} W_j (pi_j_new - pi_j)
   13.        pi_j <- pi_j_new

   14.    check convergence

   15. u_{n+1} <- u_free + du
   16. store pi_g, pi_gamma and the contact percussions for the next step

Convergence check
-----------------

The same mixed absolute/relative residual is used as for Projected Jacobi:

.. math::

   \mathrm{res}^{(k)}
   = \frac{1}{\sqrt{n_u}}
     \left\lVert
       \frac{u^{(k+1)} - u^{(k)}}
            {\varepsilon_\mathrm{abs}\,\mathbf{1} + \varepsilon_\mathrm{rel}\,\max(|u^{(k+1)}|,|u^{(k)}|)}
     \right\rVert_2 \leq 1

Nesterov acceleration and warmstart are supported with the same behaviour as
Projected Jacobi (see :doc:`projected_jacobi`).

Difference from Projected Jacobi
---------------------------------

Both solvers iterate on contact impulses and use the same preconditioner
diagonal, but:

- **PJ** updates all blocks using the velocity from the *previous* iteration
  and only then propagates the combined impulse change. This makes each sweep
  embarrassingly parallel.
- **PGS** propagates each block's impulse change to :math:`\delta u`
  immediately before moving to the next block. Later blocks therefore use more
  up-to-date information within the same sweep.

In practice PGS hopes to achieve better asymptotic runtime than solving the full constraint matrix using LU decomposition.

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
   * - ``solver.relaxation``
     - Over-relaxation factor on all solver updates
   * - ``solver.alpha``
     - Step-size multiplier for only contact impulse updates
   * - ``solver.nesterov``
     - Enable Nesterov momentum (default false)
   * - ``solver.warmstart``
     - Reuse previous-step impulses as initial guess (default true)
