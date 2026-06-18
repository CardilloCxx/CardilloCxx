Physics Pipeline
================

Every call to ``engine.step()`` drives the simulation forward by one fixed time
step. Below is a map of what happens inside that step. Each numbered box
corresponds to a section on this page.

.. contents:: Steps
   :local:
   :depth: 1

----

1 - Sync world state
--------------------

Before any math runs, the pipeline checks whether the scene has changed since
the last step. Three lazy flags are tracked on the world:

- **Structure changed** - a body or constraint was added or removed. Triggers a
  full reassignment of body indices, mass matrix rebuild, and collision scene
  rebuild.
- **State changed** - a position or velocity was set from outside the solver (e.g.
  a teleport or an initial condition). Triggers re-reading all ECS positions and
  velocities into flat global vectors.
- **Forces changed** - gravity was changed or an external force was applied.
  Triggers rebuilding the force vector, then immediately clears the
  per-entity force and torque components so they do not carry over to the next
  step.

Only the dirty paths are re-evaluated; a step with no scene changes skips all
three.

.. note::
   ``DynamicsAssembler::refreshState()`` consumes these flags and holds the flat
   numerical representation of the world: stacked position vector **q**, velocity
   vector **v**, force vector **f**, and the diagonal mass/inertia data.

----

2 - Advance positions, first half
----------------------------------

Positions and orientations are updated *explicitly* from the current velocity
using a fraction of the time step: ``(1 - theta) * dt``.

For ``theta = 1.0`` (the default, fully implicit Moreau) this term is zero and
the first half-step does nothing. For ``theta = 0.5`` (midpoint) positions move
forward by half a step before collision is evaluated.

.. note::
   The Moreau theta parameter is set via ``moreau.theta`` in the config.

----

3 - Drive trajectories and derived geometry
--------------------------------------------

Before collision sees the world, three derived-state updates run:

**Trajectories**
   Every entity carrying a trajectory component has its pose and velocity
   overwritten by evaluating the user callback at the current elapsed time.
   If only a position function was supplied, the velocity is computed by
   finite-differencing the position at ``t`` and ``t + dt``. This is the point
   where kinematic actors, conveyor slats, and scripted grippers update.

**Rigid state cache**
   A cached copy of each entity's pose is refreshed so downstream systems
   (constraints, collision) use consistent data.

**Beam segment geometry**
   Each beam segment's current length is recomputed from its neighbours'
   positions. If the length changed enough to alter the collider dimensions,
   the collision scene is flagged dirty so it will be rebuilt in the next step.

.. note::
   ``DerivedEntitySync::updateEntities`` orchestrates these three jobs;
   trajectory evaluation is in ``Trajectory::update``.

----

4 - Detect collisions
---------------------

The collision pipeline runs in three sub-stages:

**Rebuild scene (if dirty)**
   If bodies were added, removed, or beam geometry changed, the COAL broad-phase
   scene is rebuilt from scratch: collision geometry is created for every entity
   that carries ``C_Collidable``, registered with the broad-phase manager, and
   indexed for fast reverse lookup. Static geometry transforms are baked in here.

**Update transforms**
   Only dynamic entities (physics objects and trajectory-driven bodies) have
   their AABB updated each step.

**Broad phase → narrow phase → contact patches**
   The broad phase produces candidate pairs. Each pair is passed to COAL's
   narrow-phase collision check to find actual penetrating contacts and compute
   contact patches. Disabled pairs and static-static pairs are skipped.

The output is a flat list of contact structs, one per contact point, each
carrying the two entities involved, contact normal and tangent directions, body-local
attachment points, penetration depth, friction coefficient (combined from
per-body values), and the last-step impulse for warmstarting.

.. note::
   ``CollisionCoal::detectAll()`` runs the full pipeline and owns the
   authoritative contact buffer. The contact tracker inside matches current
   contacts to previous-step contacts to transfer warmstart impulses.

----

5 - Assemble contact and constraint operators
---------------------------------------------

With contacts and current poses known, the pipeline builds the numerical
operators the solver needs.

**Contact Jacobian W**
   A sparse matrix with one row per contact scalar (1 for frictionless contacts,
   3 for frictional: normal plus two tangent directions). Each row maps a body's
   velocity DOFs to the relative velocity at the contact point. Dynamic bodies
   contribute columns; static bodies contribute a constant velocity bias term
   instead.

**Constraint Jacobians Wg and W-gamma**
   Each constraint pattern (hinge, beam, distance, etc.) contributes rows that
   encode how the constraint gap changes with velocity. The split into Wg
   (spring rows) and W-gamma (damper rows) lets the solver treat position
   restoration and velocity damping separately.

**Compliance diagonals**
   Each spring row has a compliance value ``C = 1/K`` and each damper row a
   compliance ``A = 1/D``. Hard constraints (infinite stiffness) appear as zero
   compliance and are assembled as equality constraints.

**Position error**
   For compliant constraints the current position error is stored. A small
   Baumgarte correction term (scaled by ``constraint_bias_factor``) is built from
   this when forming the solver right-hand side, which gently corrects constraint
   drift over time without adding energy.

.. note::
   ``rebuildW_()`` builds the contact Jacobian. ``rebuildInteractionW_()`` builds
   constraint Jacobians and compliance data by calling ``getConstraint()`` on
   every ``ConstraintPattern``.

----

6 - Solve for new velocities
-----------------------------

The solver receives all the operators from step 5 and computes a new velocity
state that simultaneously satisfies the equations of motion, the constraint
forces, and the contact non-penetration condition.

All solvers operate on the same extended variable vector:

.. math::

   x = \begin{pmatrix} v \\ \lambda_g \\ \lambda_\gamma \end{pmatrix}

where :math:`v` is the new body velocity, :math:`\lambda_g` are spring/hard
constraint multipliers, and :math:`\lambda_\gamma` are damper multipliers.

The right-hand side encodes the current velocity, external and gyroscopic forces,
carry-over constraint multipliers from the previous step, and the Baumgarte
position correction.

Frictionless contact impulses are constrained to be non-negative (pushing
only). Frictional contacts additionally enforce the Coulomb friction cone.

Different solver backends approach the contact problem differently. See
:doc:`solvers/index` for a per-solver breakdown.

At the end of the solve, spring and damper multipliers are stored for the next
step's right-hand side, and contact impulses are stored for warmstarting.

----

7 - Write back velocities and finish positions
-----------------------------------------------

The velocity solution is scattered back into the ECS: linear velocity, angular
velocity, and derived accelerations (computed from the velocity change) are
written per entity.

Positions are then advanced through a second implicit integration, covering
the remaining ``theta * dt`` fraction of the step using the newly solved velocity.

----

8 - Write output
-----------------

After the step is fully integrated, two output writers run:

**VTK files**
   A ``.vtk`` file is written every ``output.interval_steps`` steps to the
   configured output folder. Each file contains the geometry and position of all
   entities tagged as visual. The step 0 frame is written *before* the first
   integration so the initial state is always captured.

**Tracked CSV**
   Entities tagged with ``engine.track(entity, "name")`` have their current
   position and orientation appended to a CSV file each step. This gives a
   per-step time series for selected bodies without loading the full VTK series.

Both outputs are opened lazily on first write. If ``output.interval_steps``
is 0, VTK writing is disabled entirely.
