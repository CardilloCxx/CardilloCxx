Setup and Configuration
=======================

The ``config::Config`` struct (``src/config/config.hpp``) controls every global
parameter: time step, gravity, solver, collision settings, output, and more.
A ``PhysicsEngine`` is always initialized from a ``Config``. It can be loaded via
``ConfigReader::fromFile()`` or built in code.

Each config key has a **variable name** (the ``Config`` struct field) and one or
more **file keys** that can appear in the configuration file. In config files,
key names are matched literally by ``ConfigReader``; boolean values are parsed
case-insensitively. Comments use ``#``.

.. contents:: On this page
   :local:
   :depth: 2

Loading a config
----------------

.. code-block:: cpp

   #include "config/config.hpp"
   #include "physics/api/physics.hpp"

   cardillo::config::Config cfg =
       cardillo::config::ConfigReader::fromFile("examples/scenes/pendulum/scene.config");

   cardillo::physics::PhysicsEngine engine(cfg);

You can also build ``Config`` in code and call ``engine.initFromConfig(cfg)`` on a
default-constructed engine.

All config keys
---------------

Simulation timing
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 35 12 35

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``sim_T``
     - ``sim.T``
     - ``5.0``
     - Total simulation duration in seconds. ``engine.isFinished()`` returns
       ``true`` once elapsed time reaches this value.
   * - ``sim_dt``
     - ``sim.dt``
     - ``1e-3``
     - Fixed time step in seconds (used by the no-argument ``step()``).
   * - ``sim_gravity``
     - ``sim.gravity``
     - ``(0, 0, -9.81)``
     - World-frame gravity vector (m/s²). In cfg.txt: three space- or comma-separated
       values, e.g. ``sim.gravity = 0 0 -9.81``.

Output settings
~~~~~~~~~~~~~~~

The simulation writes VTK output files that can be opened directly in
`ParaView <https://www.paraview.org/>`_ for inspection and animation.

.. list-table::
   :header-rows: 1
   :widths: 25 35 15 25

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``output_folder``
     - ``output.folder``
     - ``"./vtk_out"``
     - Output directory for VTK files.
   * - ``output_filename_prefix``
     - ``output.filename_prefix``
     - ``"scene"``
     - Filename prefix for the VTK series (e.g. ``pendulum_0001.vtp``).
   * - ``output_interval_steps``
     - ``output.interval_steps``
     - ``50``
     - Write one VTK frame every *N* time steps (minimum 1).
   * - ``output_heightfield_stride``
     - ``output.heightfield_stride``
     - ``8``
     - Decimation factor for height-field VTK triangles (minimum 1; larger reduces file size).
   * - ``output_write_contacts``
     - ``output.write_contacts``
     - ``false``
     - Also write per-step contact geometry VTK files.
   * - ``output_contacts_body_vectors``
     - ``output.contacts_body_vectors``
     - ``false``
     - Include per-contact body-frame normals/tangents/points in contact VTK.

Collision settings
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 35 15 20

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``collision_broadphase``
     - ``collision.broadphase``
     - ``"dynamic_aabb"``
     - Broad-phase algorithm. Options: ``dynamic_aabb``, ``dynamic_aabb_array``,
       ``naive``, ``sap``, ``ssap``.
   * - ``collision_disable_all``
     - ``collision.disable_all``
     - ``false``
     - Globally disable collision detection (no collidable components, no contact computation).
       Accepts ``true/1/yes/on`` (case-insensitive).
   * - ``collision_security_margin``
     - ``collision.security_margin``
     - ``0.0``
     - Extra clearance (m) added around each shape for the narrow-phase query.
       A small positive value improves robustness for thin or fast-rotating objects.
   * - ``collision_max_raw_contacts``
     - ``collision.max_raw_contacts``
     - ``1024``
     - Maximum raw contacts requested per shape pair from the narrow phase (minimum 1).
   * - ``collision_max_patches``
     - ``collision.max_patches``
     - ``4``
     - Maximum contact patches extracted per shape pair.
   * - ``collision_use_patch_vertices``
     - ``collision.use_patch_vertices``
     - ``true``
     - Emit patch-vertex contacts rather than raw contacts (more stable for planar contacts).
       Accepts ``true/1/yes/on``.
   * - ``collision_match_max_dist``
     - ``collision.match_max_dist``
     - ``0.02``
     - Maximum distance (m) between a current and previous contact for warmstarting impulse transfer.
       Clamped to ≥ 0.
   * - ``collision_min_pair_contact_distance``
     - ``collision.min_pair_contact_distance``
     - ``0.0``
     - Contacts within the same entity pair closer than this distance are deduplicated (keeping deeper penetrations).
       Set to 0 to disable; clamped to ≥ 0.

Friction settings
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 35 15 30

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``friction_enable``
     - ``friction.enable``
     - ``false``
     - Enable friction (adds tangent Jacobian rows to the contact problem).
       Accepts ``true/1/yes/on``. Roughly doubles solve cost.
   * - ``friction_combine``
     - ``friction.combine``
     - ``"min"``
     - Per-pair friction coefficient combination mode: ``min``, ``arith``, or ``geom``.
   * - ``friction_default_mu``
     - ``friction.default_mu``
     - ``0.5``
     - Default Coulomb coefficient when a body's ``RigidProps::friction`` is negative.

Debug / diagnostics
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 35 15 65

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``debug_rb``
     - ``debug.rb``
     - ``false``
     - Enable rigid-body contact / W diagnostics in the Moreau integrator.
       Accepts ``true/1/yes/on``.
   * - ``debug_pj``
     - ``debug.pj``
     - ``false``
     - Enable Projected Jacobi iteration logging. Accepts ``true/1/yes/on``.
   * - ``debug_mesh``
     - ``debug.mesh``
     - ``false``
     - Print mesh normalization info (volume, COM, inertia) at load time.
       Accepts ``true/1/yes/on``.

Solver selection and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 35 18 27

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``solver``
     - ``solver.type``
     - ``ProjectedJacobi``
     - Which contact solver to use (pj, pgs, ip, qoco, clarabel, cg)
   * - ``integrator``
     - ``integrator``, alias: ``solver.name``
     - ``Moreau``
     - The time integrator type (currently only Moreau is supported).
   * - ``pj_max_iterations``
     - ``solver.max_iterations``
     - ``200000``
     - Hard iteration cap for Projected Jacobi (minimum 1). Clamped to ≥ 1.
   * - ``pj_tol_abs``
     - ``solver.tol_abs``
     - ``1e-4``
     - Absolute convergence tolerance on the velocity residual.
   * - ``pj_tol_rel``
     - ``solver.tol_rel``
     - ``1e-4``
     - Relative convergence tolerance.
   * - ``pj_relaxation``
     - ``solver.relaxation``
     - ``0.9``
     - Over-relaxation factor ∈ (0, 1]. Below 1 stabilizes but slows convergence.
   * - ``pj_alpha``
     - ``solver.alpha``
     - ``0.3``
     - Step-size multiplier for each Jacobi update.
   * - ``pj_nesterov``
     - ``solver.nesterov``
     - ``false``
     - Enable Nesterov momentum acceleration. Accepts ``true/1/yes/on``.
   * - ``pj_nesterov_beta_threshold``
     - ``solver.nesterov_beta_threshold``
     - ``0.995``
     - Reset Nesterov momentum when beta exceeds this threshold (range 0–1).
   * - ``pj_nesterov_restart_limit``
     - ``solver.nesterov_restart_limit``
     - ``4``
     - After this many momentum resets in a solve, disable momentum for the remainder.
   * - ``pj_warmstart``
     - ``solver.warmstart``
     - ``true``
     - Reuse contact impulses from the previous step as the initial guess.
       Accepts ``true/1/yes/on``. Almost always beneficial.

Integrator (Moreau) settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 40 18 17

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``moreau_theta``
     - ``moreau.theta``, alias: ``solver.theta``
     - ``1.0``
     - Moreau theta parameter. ``1.0`` is fully implicit (most stable),
       ``0.5`` is midpoint (second-order accurate).
   * - ``moreau_implicit_gyroscopy``
     - ``moreau.implicit_gyroscopy``, alias: ``solver.implicit_gyroscopy``
     - ``false``
     - Treat the gyroscopic torque term implicitly to conserve angular momentum.
       This makes the effective mass matrix non-symmetric. Accepts
       ``true/1/yes/on``.
   * - ``moreau_lambda_theta``
     - ``moreau.lambda_theta``, alias: ``solver.lambda_theta``
     - ``false``
     - Enable θ method for integrating the Lagrange multipliers in the Moreau
       integrator. Accepts ``true/1/yes/on``.

Interior-point solver settings (QOCO / Clarabel)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters are shared across both interior-point backends and accept multiple
prefixes as shown:

.. list-table::
   :header-rows: 1
   :widths: 20 45 18 27

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``qoco_backend``
     - ``qoco.backend``
       accepts: ``auto``, ``cpu``, ``cuda``
     - ``"auto"``
     - QOCO backend selection.
   * - ``ip_kkt_static_reg``
     - ``solver.kkt_static_reg``
     - ``1e-8``
     - Static regularization applied to the KKT matrix for interior-point solvers.
   * - ``ip_kkt_dynamic_reg``
     - ``solver.kkt_dynamic_reg``
     - ``1e-8``
     - Dynamic regularization parameter for the KKT system.
   * - ``ip_iter_ref_iters``
     - ``solver.iter_ref_iters``
     - ``0``
     - Maximum number of iterative refinement passes for the KKT solve (minimum 0).

Position-error bias
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 45 18 32

   * - Variable name
     - cfg.txt key(s)
     - Default
     - Meaning
   * - ``constraint_bias_factor``
     - ``solver.constraint_bias_factor``; also ``constraint_bias_factor``,
       ``baumgarte_bias_factor``, ``baumgarte_bias``, ``baumgarte_factor``,
       ``constraint_bias``
     - ``0.001``
     - Baumgarte-style position-error correction gain. Larger values fix drift
       faster but can introduce artificial energy; 0 disables it entirely.

Scene selection
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Variable name
     - cfg.txt key(s)
     - Meaning
   * - ``scene_name``
     - ``scene.name``
     - Selects which example scene to build. If not set (default ``"none-specified"``), the engine lists available scenes and exits.

Config file syntax
------------------

Each line is ``key = value``. Lines starting with ``#`` are
comments. Inline comments after ``#`` are supported:

.. code-block:: ini

   sim.T        = 10.0          # total simulation time (seconds)
   sim.dt       = 5e-4           # time step
   sim.gravity  = 0 0 -9.81     # gx gy gz (space or comma separated)

   friction.enable    = on       # accepts true/1/yes/on, case-insensitive
   friction.combine   = geom
   friction.default_mu = 0.3

   collision.broadphase = dynamic_aabb
   collision.security_margin = 0.001
   collision.disable_all = false

   solver.max_iterations = 100000
   solver.nesterov       = true      # accepts true/1/yes/on
   solver.warmstart      = on

   solver.type        = projected_jacobi   # see table above for all aliases
   integrator         = moreau             # selects the time integrator (alias: solver.name)

   # Baumgarte bias; many cfg.txt names work interchangeably
   constraint_bias_factor = 0.01          # also works: baumgarte_bias_factor, baumgarte_bias

   # Interior-point settings accept multiple prefixes
   solver.kkt_static_reg = 1e-6           # also: ip.kkt_static_reg, qoco.kkt_static_reg
   clarabel.static_regularization_constant = 1e-6   # Clarabel-specific alias for same value

   output.folder          = ./vtk_out
   output.interval_steps  = 20
   output.filename_prefix = my_scene

   scene.name             = pendulum    # selects the example scene

The main loop
-------------

After populating the engine you drive it with a simple while-loop:

.. code-block:: cpp

   cardillo::config::Config cfg =
       cardillo::config::ConfigReader::fromFile(argv[1]);
   cardillo::physics::PhysicsEngine engine(cfg);

   myScene.populate(engine);

   real_t t  = 0.0;
   real_t dt = cfg.sim_dt;

   while (!engine.isFinished()) {
       myScene.updateScene(engine, t, dt);
       engine.step();                        // advance by cfg.sim_dt
       t += dt;
   }
