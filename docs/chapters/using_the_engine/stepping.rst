Stepping and Runtime Control
=============================

Once bodies and constraints are in place, :cpp:func:`step() <cardillo::physics::PhysicsEngine::step>` advances the simulation.
Everything else in this chapter covers runtime modifications: applying external
forces, reading and overriding state, labelling bodies for output, and freezing
or releasing bodies mid-simulation.

Advancing the simulation
------------------------

:cpp:class:`PhysicsEngine <cardillo::physics::PhysicsEngine>`
~~~~~~~~~~~~~~~~~

:cpp:func:`step() <cardillo::physics::PhysicsEngine::step>` advances by the configured time step :cpp:member:`sim_dt <cardillo::config::Config::sim_dt>`. This is the standard call inside the main loop:

.. code-block:: cpp

   while (!engine.isFinished()) {
       engine.step();
   }

``engine.step(real_t dt)``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Advances by an explicit ``dt``, ignoring ``cfg.sim_dt`` for this one step:

.. code-block:: cpp

   engine.step(5e-4);   // half a millisecond step

Useful for variable-rate drivers or when you want a single warmup step at a
different resolution before the main loop.

:cpp:func:`isFinished() <cardillo::physics::PhysicsEngine::isFinished>`
~~~~~~~~~~~~~~~~~~~~~~~

Returns ``true`` when the pipeline has executed the configured number of
steps, computed as ``ceil(cfg.sim_T / cfg.sim_dt)`` during initialization.
This matches the default ``step()`` loop. If you drive the engine with
``step(dt)`` using a different ``dt``, ``isFinished()`` still follows the
configured step count rather than accumulated wall-clock simulation time.

Gravity
-------

Gravity is set during :cpp:func:`initFromConfig() <cardillo::physics::PhysicsEngine::initFromConfig>` from ``cfg.sim_gravity``. Change it at
any time:

.. code-block:: cpp

   engine.setGravity(Vector3r(0, 0, -9.81));   // standard Earth gravity (−z)
   engine.setGravity(Vector3r::Zero());         // zero-g environment

   const Vector3r& g = engine.gravity();        // read current gravity

Applying external forces
------------------------

Force and torque are accumulated on the entity for the next assembly pass, then
zeroed by the dynamics assembler after they have been consumed. Call these from
inside ``SceneBase::updateScene`` (or before the first ``step()``) to apply
persistent loads:

.. code-block:: cpp

   // Apply an inertial-frame force and torque to entity e
   engine.applyForce(e,
       Vector3r(0, 0, 100),    // force in N (inertial frame)
       Vector3r::Zero()        // torque in Nm (inertial frame)
   );

   // Apply a pure torque defined in inertial frame
   engine.applyInertialTorque(e, Vector3r(0, 0, 5));

.. note::
   ``applyForce`` and ``applyInertialTorque`` **accumulate**; calling them
   multiple times per step adds up. Forces are cleared automatically after
   each ``step()`` call by the integrator.

Reading state
-------------

.. code-block:: cpp

   VectorXr q = engine.getPosition(e);        // generalised position (7 for RB)
      MatrixXXr M = engine.getMass(e);            // 6×6 or 3×3 mass matrix
      Vector3r  I = engine.getInertiaDiag(e);     // Body Basis inertia diagonal
      real_t   Ek = engine.getKineticEnergy(e);   // translational + rotational KE

   For rigid bodies the generalised position vector is
``[px, py, pz, qx, qy, qz, qw]`` because it is returned from
:cpp:type:`Quaternion4r <cardillo::Quaternion4r>`.

Setting state
-------------

Use these to teleport a body or impose an initial condition mid-simulation.
Each setter marks the world state dirty so the pipeline re-assembles on the
next ``step()``:

.. code-block:: cpp

   engine.setPosition(e, Vector3r(1, 0, 2));
   engine.setOrientation(e, Quaternion4r::Identity());
   engine.setLinearVelocity(e, Vector3r(0, 0, -1));
   engine.setAngularVelocity(e, Vector3r(0, 1, 0));

   // Set both linear and angular velocity by applying an impulse (force over one time step)
   engine.setVelocityByForce(e, Vector3r(0, 0, -1), Vector3r(0, 1, 0));

.. warning::
   Teleporting a body that is connected to other bodies by hard constraints
   can introduce large position errors. Either mark the structure dirty and
   let the Baumgarte correction handle it, or rebuild constraints after the
   move.

Marking the world dirty
-----------------------

The pipeline uses three lazy-rebuild flags that are set automatically by most
mutation methods. You can also set them manually:

.. code-block:: cpp

   engine.markStructureDirty();  // body added/removed; rebuild mass matrix etc.

This is rarely needed manually because every body-creation and constraint method
calls it internally.

Tracking bodies for output
--------------------------

Attach a named label to any entity to include its state over time in a CSV file for plotting.

.. code-block:: cpp

   engine.track(arm, "pendulum_arm");
   engine.track(tip, "beam_tip");

Freezing a body mid-simulation
-------------------------------

Convert a dynamic body into a static obstacle without removing it:

.. code-block:: cpp

   engine.makeStatic(body);

This removes the entity's dynamic physics components (including
:cpp:struct:`C_PhysicsObject <cardillo::C_PhysicsObject>`, rigid/point-mass tags, body index, mass, inertia, and any
queued external force/torque data) and marks the structure dirty. Useful for
staged simulations where a body should be settled, locked in place, and then
other bodies dropped on top.

Accessing the ECS registry
---------------------------

For custom queries or one-off component reads that are not exposed by the
:cpp:class:`PhysicsEngine <cardillo::physics::PhysicsEngine>` API, access the underlying EnTT registry directly:

.. code-block:: cpp

   const entt::registry& reg = engine.ecs();

   // Example: iterate all bodies that have a friction component
   reg.view<cardillo::C_Friction>().each([](auto e, const auto& f) {
       std::cout << "mu = " << f.mu << "\n";
   });

.. caution::
   Modifying ECS components directly (adding, removing, or overwriting them)
   bypasses the engine's dirty flags. If you change any component that affects
   assembly, call :cpp:func:`markStructureDirty() <cardillo::physics::PhysicsEngine::markStructureDirty>` afterwards.
