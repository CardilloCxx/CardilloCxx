Time-Stepping for Non-Smooth Mechanical Systems
===============================================

This chapter details the discrete equations used in Cardillo's Moreau-style 
time-stepping integrator for rigid bodies subject to unilateral constraints and 
friction :cite:`moreau1988`. The integrator employs a **velocity-level splitting** scheme with a generalized :math:`\theta`-method for velocity updates.

Non-smooth mechanical systems are characterized by abrupt velocity changes and impulsive reactions. The Moreau time-stepping scheme integrates the 
system at the velocity level, which naturally accommodates impulsive constraints while maintaining stability during smooth phases.

Symbols and Notation
--------------------

* :math:`\mathbf{M}`: Block-diagonal mass/inertia matrix.
* :math:`\mathbf{u}`: Stacked generalized velocity vector.
* :math:`\mathbf{W}_g, \mathbf{W}_\gamma`: Constraint Jacobian blocks (spring and damper rows).
* :math:`\mathbf{W}_{NT}`: Contact Jacobian (stacked normal and tangential rows).
* :math:`\mathbf{\Lambda}_{NT}`: Accumulated contact impulses over the step.
* :math:`\boldsymbol{\gamma}_{NT}^+`: Post-impact contact velocity in Normal and Tangential directions.
* :math:`\mathbf{f}^{\mathrm{ext}}`: External forces (e.g. gravity).
* :math:`\mathbf{v}_g, \mathbf{v}_\gamma`: Velocity-source terms from trajectory-driven bodies.
* :math:`\mathbf{C}, \mathbf{A}`: Compliance/attenuation diagonals.
* :math:`\mathbf{g}`: Position-level measures of the spring rows (constraint violation, rod strains).
* :math:`\boldsymbol\lambda_g, \boldsymbol\lambda_\gamma`: Lagrange multiplier densities.
* :math:`\mathbf{\Lambda}_g, \mathbf{\Lambda}_\gamma`: Solver-scaled impulses, defined as :math:`\mathbf{\Lambda}_g = -\Delta t\,\boldsymbol\lambda_{g,\theta}` and :math:`\mathbf{\Lambda}_\gamma = -\Delta t\,\boldsymbol\lambda_{\gamma,\theta}`.
* :math:`\theta, \Delta t`: Integration parameter and time-step size.

Continuous Equations of Motion
------------------------------

The system is governed by the following measure-differential equations:

.. math::
    \begin{aligned}
      \dot{\mathbf{q}} &= \mathbf{B}(\mathbf{q}) \mathbf{u} \\
      \mathbf{M} \operatorname{d}\!{\mathbf{u}} &= \bigl[\mathbf{f}^\mathrm{ext}(t, \mathbf{q}) + \mathbf{f}^\mathrm{gyr}(\mathbf{u})\bigr] \operatorname{d}\!t + \mathbf{W}_g(\mathbf{q}) \operatorname{d}\!{\boldsymbol{\pi}_g} + \mathbf{W}_\gamma(\mathbf{q}) \operatorname{d}\!{\boldsymbol{\pi}_\gamma} + \mathbf{W}_N(\mathbf{q}) \operatorname{d}\!\boldsymbol{\pi}_N + \mathbf{W}_T(\mathbf{q}) \operatorname{d}\!\boldsymbol{\pi}_T \\
      \mathbf{0} &= \mathbf{C} \boldsymbol{\lambda}_g + \mathbf{g}(t, \mathbf{q}) \\
      \mathbf{0} &= \mathbf{A} \boldsymbol{\lambda}_\gamma + \mathbf{W}_\gamma^\top(\mathbf{q}) \mathbf{u}
    \end{aligned}

Velocity-Level Splitting
------------------------

The time step is decomposed into three phases using an intermediate configuration 
:math:`\mathbf{q}_{n+\theta}` at the intermediate time :math:`t_{n+\theta} = t_n + (1-\theta)\,\Delta t`.
Kinematically driven bodies (trajectories) and prescribed rates of constraints
(e.g. a driven rest length) are advanced to the same time before the measures
:math:`\mathbf{g}` and Jacobians :math:`\mathbf{W}` are evaluated, and by the
remaining :math:`\theta\,\Delta t` after the second drift.

**1. First Phase (Explicit Position Update):**
Updates configuration based on the current velocity:

.. math::
    \mathbf{q}_{n+\theta} = \mathbf{q}_n + (1-\theta)\,\Delta t\,\mathbf{B}(\mathbf{q}_n)\,\mathbf{u}_n

**2. Second Phase (Implicit Velocity Update) at** :math:`n+\theta`:
The momentum balance is solved using a generalized :math:`\theta`-method. 

.. note::
    While :math:`\theta=0.5` is often desired for second-order accuracy, 
    choosing :math:`\theta \in (0.5, 1]` introduces numerical damping, which 
    improves stability for stiff mechanical systems (e.g., beam elements).

Approximating the gyroscopic term as :math:`\mathbf{f}^{\mathrm{gyr}} \approx \mathbf{G}(\mathbf{u}_n)\mathbf{u}_{n+1}` and defining 
:math:`\mathbf{M}_{\mathrm{eff}} := \mathbf{M} - \Delta t\,\mathbf{G}(\mathbf{u}_n)`, we use the backward Euler 
approximation for impulses:

.. math::
   \Delta \boldsymbol{\pi}_{g,n+1} = -\Delta t\,\boldsymbol{\lambda}_{g,n+1}, \quad
   \Delta \boldsymbol{\pi}_{\gamma,n+1} = -\Delta t\,\boldsymbol{\lambda}_{\gamma,n+1}.

The discrete momentum equation becomes:

.. math::
    \mathbf{M}_{\mathrm{eff}}\mathbf{u}_{n+1} = \mathbf{M}\mathbf{u}_n + \Delta t\,\mathbf{f}^{\mathrm{ext}}_{n+\theta} + \mathbf{W}_{NT}\,\boldsymbol{\Lambda}_{NT} + \mathbf{W}_g \boldsymbol{\Lambda}_g + \mathbf{W}_\gamma \boldsymbol{\Lambda}_\gamma

**3. Third Phase (Second Drift):**
Updates configuration for the next time step:

.. math::
    \mathbf{q}_{n+1} = \mathbf{q}_{n+\theta} + \theta\,\Delta t\,\mathbf{B}(\mathbf{q}_{n+\theta})\,\mathbf{u}_{n+1}

Discrete Constraint and Contact Laws
------------------------------------

The unknowns of the force laws are the :math:`\theta`-weighted forces
:math:`\boldsymbol\lambda_{g,\theta}` and :math:`\boldsymbol\lambda_{\gamma,\theta}`. The spring law is
enforced on position level, linearized at the intermediate configuration
(the only force law of the implementation):

.. math::
    \mathbf{C}\,\boldsymbol{\lambda}_{g,\theta} + \mathbf{g}(t_{n+\theta}, \mathbf{q}_{n+\theta}) + \Delta t\,\bigl(\theta^2\,\mathbf{W}_{g}^\top \mathbf{u}_{n+1} - (1-\theta)^2\,\mathbf{W}_{g}^\top \mathbf{u}_n\bigr) + (2\theta - 1)\,\Delta t\,\mathbf{v}_g = \mathbf{0}

It is the first-order expansion about :math:`(t_{n+\theta}, \mathbf{q}_{n+\theta})` of
:math:`\mathbf{C}\boldsymbol\lambda_{g,\theta} + \mathbf{g}(t_n, \mathbf{q}_n) + \theta\,\Delta t\,\bigl(\mathbf{W}_g^\top \mathbf{u}_{n+\theta} + \mathbf{v}_g\bigr) = \mathbf{0}`
with :math:`\mathbf{u}_{n+\theta} = (1-\theta)\,\mathbf{u}_n + \theta\,\mathbf{u}_{n+1}`, which
would require an additional evaluation of :math:`\mathbf{g}` at the beginning of
the step. No multiplier of
the previous step appears, so no initial multipliers are needed, neither for
compliant elements nor for perfect constraints (:math:`\mathbf{C} = \mathbf{0}`), and the
constraints do not drift. Perfect constraints require :math:`\theta > 0.5`. The
damper law is

.. math::
    \mathbf{A}\,\boldsymbol{\lambda}_{\gamma,\theta} + \mathbf{W}_{\gamma}^\top\bigl((1-\theta)\,\mathbf{u}_n + \theta\,\mathbf{u}_{n+1}\bigr) + \mathbf{v}_\gamma = \mathbf{0}

For contact constraints with friction, the following complementarity/feasibility conditions apply at each contact :math:`j`:

.. math::
    \begin{aligned}
        &0 \le \gamma_{N,j}^+ \ \perp\ \Lambda_{N,j} \ge 0, \\
        &\lVert\boldsymbol{\Lambda}_{T,j}\rVert \le \mu_j\,\Lambda_{N,j}, \\
        &\boldsymbol{\gamma}_{T,j}=\mathbf{0}\ \Rightarrow\ \lVert\boldsymbol{\Lambda}_{T,j}\rVert \le \mu_j\,\Lambda_{N,j}\qquad\text{(sticking)}, \\
        &\boldsymbol{\gamma}^+_{T,j}\neq\mathbf{0}\ \Rightarrow\ \boldsymbol{\Lambda}_{T,j} = -\mu_j\,\Lambda_{N,j}\,\dfrac{\boldsymbol{\gamma}^+_{T,j}}{\lVert\boldsymbol{\gamma}^+_{T,j}\rVert} \quad\text{and}\quad \lVert\boldsymbol{\Lambda}_{T,j}\rVert = \mu_j\,\Lambda_{N,j}\qquad\text{(sliding)}.
    \end{aligned}

Scaled Block Linear System
--------------------------

With the impulses :math:`\mathbf{\Lambda}_g = -\Delta t\,\boldsymbol\lambda_{g,\theta}` and :math:`\mathbf{\Lambda}_\gamma = -\Delta t\,\boldsymbol\lambda_{\gamma,\theta}`,
the discrete equations can be rearranged into the following scaled block linear system:

.. math::
    \underbrace{%
    \begin{pmatrix} 
    \mathbf{M}_{\mathrm{eff}} & \mathbf{W}_g & \mathbf{W}_{\gamma} \\
    \mathbf{W}_g^\top & -\dfrac{1}{\theta^2\,\Delta t^2}\mathbf{C} & \mathbf{0} \\
    \mathbf{W}_{\gamma}^\top & \mathbf{0} & -\dfrac{1}{\theta\,\Delta t}\mathbf{A} 
    \end{pmatrix}%
    }_{\mathcal{S}(\mathbf{q}_{n+\theta},\Delta t,\theta)}\;
    \underbrace{%
    \begin{pmatrix} 
    \mathbf{u}_{n+1} \\
    \boldsymbol{\Lambda}_{g} \\
    \boldsymbol{\Lambda}_{\gamma} 
    \end{pmatrix}%
    }_{\mathbf{x}_{n+1}} 
    =
    \underbrace{%
    \begin{pmatrix} 
    \mathbf{M}\mathbf{u}_n + \Delta t\,\mathbf{f}^{\mathrm{ext}}_{n+\theta} + \mathbf{W}_{NT}\boldsymbol{\Lambda}_{NT} \\
    -\dfrac{1}{\theta^2\,\Delta t}\mathbf{g}(t_{n+\theta}, \mathbf{q}_{n+\theta}) + \Bigl(\dfrac{1-\theta}{\theta}\Bigr)^2\mathbf{W}_g^\top \mathbf{u}_n - \dfrac{2\theta-1}{\theta^2}\mathbf{v}_g \\
    -\dfrac{1-\theta}{\theta}\mathbf{W}_\gamma^\top \mathbf{u}_n - \dfrac{1}{\theta}\mathbf{v}_\gamma 
    \end{pmatrix}%
    }_{\mathbf{b}_n(t_{n+\theta}, \mathbf{q}_{n+\theta},\mathbf{u}_n,\Delta t,\theta,\boldsymbol{\Lambda}_{NT})}

The right-hand side contains no multipliers of the previous step; the stored
impulses only serve as warm start of the iterative solvers. The different
solver backends all solve this same linear system in some form; the
compliance factors and the spring and damper biases are computed once in
``DynamicsAssembler::springComplianceScale()``, ``springBias()`` and
``damperBias()``.

References
----------

.. bibliography:: ../refs.bib 
    :style: unsrt