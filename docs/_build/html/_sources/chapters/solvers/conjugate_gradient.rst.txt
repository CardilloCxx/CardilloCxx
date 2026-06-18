Conjugate Gradient
==================

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

The Conjugate Gradient solver applies a preconditioned CG iteration to the
Schur complement system that arises from the Moreau time step. It is
constraint-only: if any contacts are active at call time the solver throws a
runtime error. Use it for scenes that rely entirely on joints, beams, and
springs with no dynamic collisions.

The Schur complement system
----------------------------

Starting from the same extended Moreau system as PGS (see
:doc:`projected_gauss_seidel`), eliminating :math:`u_{n+1}` by block
elimination gives:

.. math::

   \mathbf{S}\,\hat\Lambda = \tilde{b},
   \qquad
   \mathbf{S} = \hat{W}^\top M^{-1} \hat{W} - \hat{C}

where:

- :math:`\hat\Lambda = (\Lambda_g^\top, \Lambda_\gamma^\top)^\top` are the
  stacked constraint multipliers
- :math:`\hat{W} = (W_g, W_\gamma)` is the stacked constraint Jacobian
- :math:`\hat{C} = \mathrm{diag}(-C/(\theta\,dt^2),\; -A/(\theta\,dt))` is
  the (negative) block-diagonal compliance matrix
- :math:`\tilde{b}` is the reduced right-hand side from the Schur complement
  of the full system

Because :math:`M^{-1}` is symmetric positive definite and the compliance
diagonals :math:`-C/(\theta\,dt^2)` and :math:`-A/(\theta\,dt)` are
non-positive, the full Schur complement :math:`\mathbf{S}` is symmetric
positive semi-definite. This makes it amenable to the conjugate gradient
method.

For hard constraints (:math:`C \to 0`, :math:`A \to 0`), the compliance
diagonal vanishes and :math:`\mathbf{S}` reduces to the classical Delassus
operator :math:`\hat{W}^\top M^{-1} \hat{W}`.

The preconditioner
------------------

The block-diagonal part of :math:`\mathbf{S}` is used as a preconditioner.
For each constraint pattern :math:`j`, the preconditioner block is:

.. math::

   D_j = \hat{W}_j^\top M^{-1} \hat{W}_j - \hat{C}_j

This is the same per-constraint block diagonal assembled by the PGS solver.
Preconditioning with :math:`D^{-1}` clusters the eigenvalues of
:math:`D^{-1}\mathbf{S}` near 1, which is the condition for CG to converge
quickly.

The matrix-vector product :math:`\mathbf{S}\,p` is evaluated implicitly
without ever forming :math:`\mathbf{S}`:

.. math::

   u_\mathrm{tmp} = M^{-1} \hat{W}\,p,
   \qquad
   \mathbf{S}\,p = \hat{W}^\top u_\mathrm{tmp} - \hat{C}\,p

This requires one sparse matrix-vector multiply with :math:`\hat{W}` and one
with :math:`\hat{W}^\top` per CG iteration.

The PCG algorithm
-----------------

.. code-block:: text

   Algorithm: Preconditioned Conjugate Gradient for Contact Impulses
   ----------------------------------------------------------------
   Input : ˜b, Ŵ, Ĉ, D⁻¹, M⁻¹, λ_prev, k_max, tol_abs
   Output: u_{n+1}

   1.  λ ← λ_prev
   2.  u_corr ← M⁻¹ Ŵ λ
   3.  r ← ˜b − (Ŵᵀ u_corr − Ĉ λ)
   4.  z ← D⁻¹ r;   p ← z;   rz ← rᵀ z

   5.  for k = 0 to k_max do
   6.      u_tmp ← M⁻¹ Ŵ p
   7.      Sp    ← Ŵᵀ u_tmp − Ĉ p
   8.      α     ← rz / (pᵀ Sp)

   9.      λ ← λ + α p
   10.     r ← r − α Sp

   11.     if ||r|| < tol_abs then break

   12.     z ← D⁻¹ r
   13.     rz_new ← rᵀ z
   14.     β ← rz_new / rz
   15.     p ← z + β p
   16.     rz ← rz_new

   17. u_corr ← M⁻¹ Ŵ λ
   18. u_{n+1} ← u_free − u_corr
   19. Store (λ, λ_prev) for next time-step

Each iteration costs two sparse matrix-vector products (with :math:`\hat{W}`
and :math:`\hat{W}^\top`) plus two preconditioner applications. For large
sparse constraint systems this is much cheaper per iteration than a direct
method. Even though the iteration should converge in at most :math:`n_\lambda` steps (the number
of multiplier unknowns), for stiff problems even with a preconditioner it might take longer.

Warmstart
---------

The multipliers from the previous step are used to initialise
:math:`\hat\Lambda^{(0)}`, which makes the initial residual small when the
constraint forces change slowly between steps. This is the primary source of
fast convergence in quasi-static or slowly varying scenes.

Limitations
-----------

- No contact support. The solver raises a runtime error if any contacts are
  present. Use PJ or PGS for scenes with collisions.
- Only symmetric problems. Implicit gyroscopy (``moreau.implicit_gyroscopy``)
  breaks the symmetry of :math:`M_\mathrm{eff}` and is not supported; the
  solver falls back to treating the system as symmetric.
- Works best for scenes where the constraint network is large but well
  structured (beam assemblies, cable networks, kinematic chains) and contacts
  are disabled.
