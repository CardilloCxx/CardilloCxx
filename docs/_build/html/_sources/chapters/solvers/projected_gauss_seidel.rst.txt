Projected Gauss-Seidel
======================

.. contents:: On this page
   :local:
   :depth: 2

The core problem
----------------

The Moreau time-stepping scheme produces the following scaled linear system at
each step (see :doc:`../physics_pipeline`, step 6):

.. math::

   \underbrace{\begin{pmatrix}
     M_\mathrm{eff} & W_g & W_\gamma \\
     W_g^\top & -\tfrac{C}{\theta\,dt^2} & 0 \\
     W_\gamma^\top & 0 & -\tfrac{A}{\theta\,dt}
   \end{pmatrix}}_{\mathcal{S}}
   \begin{pmatrix} u_{n+1} \\ \Lambda_g \\ \Lambda_\gamma \end{pmatrix}
   =
   \begin{pmatrix}
     M u_n + dt\, f^\mathrm{ext} + W_{NT}\,\Lambda_{NT} \\
     -\tfrac{C\,\Lambda_{g,n}}{\theta\,dt^2} - \tfrac{1-\theta}{\theta}W_g^\top u_n - \tfrac{1}{\theta}v_g \\
     -\tfrac{1-\theta}{\theta}W_\gamma^\top u_n - \tfrac{1}{\theta}v_\gamma
   \end{pmatrix}

where:

- :math:`u` is the stacked generalized velocity
- :math:`\Lambda_g = -dt\,\lambda_g` are scaled spring multipliers
- :math:`\Lambda_\gamma = -dt\,\lambda_\gamma` are scaled damper multipliers
- :math:`M_\mathrm{eff} = M - dt\,G(u_n)` is the effective mass including the optional implicit gyroscopic correction
- :math:`C = K^{-1}` and :math:`A = D^{-1}` are the compliance and attenuation diagonals
- :math:`v_g`, :math:`v_\gamma` are velocity source terms from trajectory-driven bodies
- :math:`\Lambda_{NT}` is the stacked contact impulse vector with normal and tangential components

Reducing to the constraint multiplier space
-------------------------------------------

Eliminating :math:`u_{n+1}` from the second and third block rows via the Schur
complement gives a reduced system purely in the constraint multipliers
:math:`\hat\Lambda = (\Lambda_g^\top, \Lambda_\gamma^\top)^\top`:

.. math::

   \mathbf{S}\,\hat\Lambda = \tilde{b},

where the Schur complement matrix is:

.. math::

   \mathbf{S} = \hat{W}^\top M_\mathrm{eff}^{-1} \hat{W} - \hat{C},

with :math:`\hat{W} = (W_g, W_\gamma)` and block-diagonal
:math:`\hat{C} = \mathrm{diag}(-C/(\theta\,dt^2),\; -A/(\theta\,dt))`.
The right-hand side is:

.. math::

   \tilde{b} = \hat{W}^\top M_\mathrm{eff}^{-1} b_u - \hat{b}

where :math:`b_u` and :math:`\hat{b}` are the corresponding blocks of the full
right-hand side. Once :math:`\hat\Lambda` is known, the velocity is recovered as:

.. math::

   u_{n+1} = \underbrace{M_\mathrm{eff}^{-1} b_u}_{u_\mathrm{free}}
           - \underbrace{M_\mathrm{eff}^{-1} \hat{W}\,\hat\Lambda}_{u_\mathrm{corr}}

The matrix-vector product :math:`\mathbf{S}\hat\Lambda` is evaluated through
the velocity correction:

.. math::

   u_\mathrm{corr} = M_\mathrm{eff}^{-1}\hat{W}\hat\Lambda,
   \qquad
   \mathbf{S}\hat\Lambda = \hat{W}^\top u_\mathrm{corr} - \hat{C}\hat\Lambda.

This avoids ever forming :math:`\mathbf{S}` explicitly.

Extending to contacts
---------------------

Adding contacts appends the contact Jacobian :math:`W_{NT}` as additional rows.
The extended system becomes:

.. math::

   \begin{pmatrix}
     M_\mathrm{eff} & \hat{W} & W_{NT} \\
     \hat{W}^\top & -\hat{C} & 0 \\
     W_{NT}^\top & 0 & 0
   \end{pmatrix}
   \begin{pmatrix} u_{n+1} \\ \hat\Lambda \\ \Lambda_{NT} \end{pmatrix}
   =
   \begin{pmatrix} b_u \\ \hat{b} \\ b_{NT} \end{pmatrix}

The zero block in the lower right means contact rows have no compliance term.
Contact impulses are not determined by this linear system alone; they are subject
to the non-convex contact projection (Signorini + Coulomb):

- **frictionless**: :math:`\Lambda_{N} \geq 0`
- **frictional**: :math:`\|\Lambda_{T}\| \leq \mu\,\Lambda_{N}`

The full multiplier vector in the implementation is therefore:

.. math::

   \lambda = \begin{pmatrix} \Lambda_g \\ \Lambda_\gamma \\ \Lambda_{NT} \end{pmatrix}

The block diagonal preconditioner
----------------------------------

PGS preconditions each block :math:`j` with the diagonal of the Schur complement
restricted to that block:

.. math::

   D_j = \hat{W}_j^\top M_\mathrm{eff}^{-1} \hat{W}_j - \hat{C}_j

For constraint patterns (hinges, beams, etc.) this is a small dense matrix
(up to 6x6) that may be inverted implicitly using a Cholesky factorization.

The PGS algorithm
-----------------

PGS differs from Projected Jacobi in that it updates each block sequentially
and propagates the velocity correction :math:`u_\mathrm{corr}` immediately,
so later blocks in the same sweep already see the effect of earlier updates.

.. code-block:: text

   Algorithm: Projected Gauss-Seidel for Contact Dynamics
   -------------------------------------------------------
   Input : rhs, W, M^{-1}, D^{-1}, k_max, relaxation
   Output: u_{n+1}

   1. Initialize λ ← λ_warm_start
   2. u_corr ← M^{-1} W^T λ

   3. for k = 0 to k_max do
   4.     for each block j in order do
   5.         r_j ← rhs_j - W_j u_corr - C_j λ_j
   6.         Δλ_j ← ω · D_j^{-1} · r_j
   7.         if j is contact block then Δλ_j ← α · Δλ_j

   8.         λ_j ← λ_j + Δλ_j

   9.         if j is frictionless then
   10.            λ_j ← min(λ_j, 0)
   11.        else if j is frictional then
   12.            λ_j ← ProjectToCoulombDisk(λ_j)

   13.        Δλ_j_proj ← λ_j - (λ_j_old)
   14.        u_corr ← u_corr + M^{-1} W_j^T Δλ_j_proj

   15.    check convergence (||Δλ|| < ε)

   16. u_{n+1} ← u_free - u_corr
   17. store λ and warm-start data for next step


.. note::
    The implementation stores normal contact impulses with
    a negative sign convention inherited from the way the contact Jacobian is
    assembled.

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
- **PGS** propagates each block's impulse change to :math:`u_\mathrm{corr}`
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
