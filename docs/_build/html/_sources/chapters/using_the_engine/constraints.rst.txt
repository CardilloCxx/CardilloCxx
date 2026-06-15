Constraints and Joints
======================

Constraints couple two bodies (or a body to the world) by restricting
translation and/or rotation. Every factory method on ``PhysicsEngine``
delegates to ``ConstraintFactory`` (``src/physics/constraints/``), inserts
a ``ConstraintPattern`` into the world's constraint list, and marks the
world structure dirty so the next ``step()`` rebuilds the solver matrices.

All factory methods return a ``size_t`` constraint index. Keep this index if
you need to drive the constraint velocity at runtime with
``setConstraintScalarVelocity`` / ``setConstraintLinearVelocity`` /
``setConstraintAngularVelocity``.

.. contents:: On this page
   :local:
   :depth: 2

JointFrame (specifying where and how a joint attaches)
-------------------------------------------------------

``JointFrame`` (``src/physics/constraints/constraints.hpp``) describes the
position and orientation of the joint in some reference frame:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Field
     - Default
     - Meaning
   * - ``r_refJ``
     - ``(0, 0, 0)``
     - Joint origin expressed in the reference frame (or in world frame when
       ``ref`` is empty).
   * - ``A_refJ``
     - ``Identity``
     - 3×3 rotation matrix defining the joint's local axes in the reference
       frame. **Column 0 is the primary (hinge) axis** used by revolute joints.
       Columns 1 and 2 are the blocked rotation axes.
   * - ``ref``
     - ``nullopt``
     - Optional entity whose current frame defines ``r_refJ`` and ``A_refJ``.
       When empty the values are taken directly in world coordinates.

The most convenient constructor is ``JointFrame::fromAxis``:

.. code-block:: cpp

   // Build a joint frame at 'pivot' whose x-column equals 'axis'.
   // Columns 1 and 2 are filled in automatically as an orthonormal basis.
   JointFrame jf = JointFrame::fromAxis(Vector3r pivot, Vector3r axis,
                                        std::optional<entt::entity> ref = {});

   // Hinge at world origin, rotating around the Y axis
   JointFrame jf = JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitY());

   // Hinge at a world point, in body A's local frame
   JointFrame jf = JointFrame::fromAxis(
       Vector3r(0.02, 0.0, 0.035),   // pivot in body A's frame
       Vector3r::UnitX(),            // hinge axis in body A's frame
       bodyA                         // reference entity
   );

You can also build a ``JointFrame`` directly from a world position and
rotation matrix:

.. code-block:: cpp

   JointFrame jf(Vector3r pivot, Matrix33r orientation);

.. note::
   When ``ref`` is set to an entity, ``r_refJ`` and the columns of ``A_refJ``
   are expressed in that entity's **current** frame at constraint-creation time.
   Changes to the entity's pose after creation do not retroactively update the
   stored attachment offsets.

How stiffness and damping work
-------------------------------

All constraints except the beam use a **compliant constraint** formulation.
Stiffness ``K`` and damping ``D`` are per-axis vectors in the joint frame.

- ``K = infinity`` (the default) → **hard** constraint: the DOF is fully
  locked at assembly time and enforced with an impulsive force.
- ``K = 0`` → **free**: no restoring force, the DOF is unconstrained.
- A finite positive ``K`` → **spring**: a restoring force proportional to the
  displacement from the rest position, and a damping force proportional to the
  relative velocity.

The compliance ``C = 1/K`` (and analogously for damping ``A = 1/D``) enters the
Delassus operator directly; no special casing needed.

Constraint types
----------------

addLinearDistanceConstraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maintains the Euclidean distance between two attachment points (one on each
body) at its value at construction time. This is a single scalar constraint
row (1 DOF).

.. code-block:: cpp

   size_t idx = engine.addLinearDistanceConstraint(
       entt::entity a,                      // body A (or entt::null for world)
       entt::entity b,                      // body B
       Vector3r rA_local = Vector3r::Zero(), // attachment offset in body A frame
       Vector3r rB_local = Vector3r::Zero(), // attachment offset in body B frame
       real_t stiffness  = infinity,         // spring stiffness (N/m); inf = hard
       real_t damping    = 0.0              // damper coefficient (Ns/m)
   );

   // Hard cable between body centres
   engine.addLinearDistanceConstraint(a, b);

   // Soft spring at offset attachment points
   engine.addLinearDistanceConstraint(a, b,
       Vector3r(0, 0, 0.5),   // attach at top of A
       Vector3r(0, 0, -0.5),  // attach at bottom of B
       1e4,                   // 10 kN/m
       50.0                   // 50 Ns/m damping
   );

The rest length ``L₀`` is computed at construction from the initial distance
between the two attachment points. Drive the constraint with a constant
scalar velocity using ``engine.setConstraintScalarVelocity(idx, v)``.

addRigidConstraint
~~~~~~~~~~~~~~~~~~

Fully welds body ``b`` to body ``a`` (or to the world when ``b`` is
``entt::null``). Equivalent to locking all 6 DOF (3 translation + 3 rotation)
with infinite stiffness. The joint position is taken as body ``a``'s current
world position at call time.

.. code-block:: cpp

   size_t idx = engine.addRigidConstraint(
       entt::entity a,             // reference body (or entt::null for world)
       entt::entity b = entt::null // body to pin (null → pins 'a' to world)
   );

   engine.addRigidConstraint(wall, bracket); // weld bracket to wall

addHingeConstraint
~~~~~~~~~~~~~~~~~~

The most commonly used joint. Allows rotation around a single axis (the
x-column of ``frame.A_refJ``) while locking all other DOF. An optional torsional
spring/damper can be added around the hinge axis.

.. code-block:: cpp

   size_t idx = engine.addHingeConstraint(
       entt::entity a,                              // first body
       entt::entity b,                              // second body (or entt::null)
       const JointFrame& frame,                     // pivot + hinge axis
       real_t K_axis  = 0.0,                        // torsional spring (Nm/rad)
       real_t D_axis  = 0.0,                        // torsional damping (Nms/rad)
       Vector3r K_trans = Vector3r::Constant(inf),  // translational stiffness
       Vector3r D_trans = Vector3r::Zero()          // translational damping
   );

Parameters:

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Parameter
     - Default
     - Meaning
   * - ``a``, ``b``
     - auto
     - The two bodies. Pass ``entt::null`` for ``b`` to pin body ``a`` to the
       world (useful for rollers, wall-mounted hinges, etc.).
   * - ``frame``
     - required
     - Defines the pivot world position and the hinge axis (col 0 of
       ``A_refJ``). Build it with ``JointFrame::fromAxis(pivot, axis)``.
   * - ``K_axis``
     - ``0``
     - Torsional spring stiffness around the hinge axis (Nm/rad). ``0`` = free
       hinge; use a large value (e.g. 1e4) for a stiff rotational spring.
   * - ``D_axis``
     - ``0``
     - Torsional damping around the hinge axis (Nms/rad). Light damping (e.g.
       0.1–1.0) prevents energy build-up in long simulations.
   * - ``K_trans``
     - ``inf`` × 3
     - Translational spring stiffness in all three joint-frame directions (N/m).
       Default infinite = rigid; lower values create a compliant pivot.
   * - ``D_trans``
     - ``(0, 0, 0)``
     - Translational damping (Ns/m).

.. note::
   Internally ``addHingeConstraint`` creates a ``TranslationRotationConstraint``
   with ``K_rot = (K_axis, inf, inf)``. The infinite stiffness on the two
   transverse rotation axes is what locks out all rotations except the hinge
   axis.

.. code-block:: cpp

   // Free hinge at world origin, rotating around Y
   JointFrame frame = JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitY());
   engine.addHingeConstraint(base, arm, frame);

   // Spring hinge with damping; mimics a torsional rubber bushing
   engine.addHingeConstraint(base, arm, frame, 200.0, 5.0);

   // Wall-mounted roller (body pinned to world)
   engine.addHingeConstraint(roller, entt::null,
       JointFrame::fromAxis(rollerCenter, Vector3r::UnitY()),
       0.0, 0.001);  // near-frictionless rotation

addTranslationalConstraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Locks all translational DOF while leaving all three rotational DOF free
(a ball-in-socket joint with translation removed).

.. code-block:: cpp

   size_t idx = engine.addTranslationalConstraint(
       entt::entity a,
       entt::entity b,
       const JointFrame& frame,
       Vector3r K_trans = Vector3r::Constant(infinity),
       Vector3r D_trans = Vector3r::Zero()
   );

This is a ``TranslationRotationConstraint`` with ``K_rot = (0, 0, 0)`` and
``K_trans`` set to the given values.

addRotationConstraint
~~~~~~~~~~~~~~~~~~~~~

Locks all rotational DOF while leaving all three translational DOF free.

.. code-block:: cpp

   size_t idx = engine.addRotationConstraint(
       entt::entity a,
       entt::entity b,
       const JointFrame& frame,
       Vector3r K_rot = Vector3r::Constant(infinity),
       Vector3r D_rot = Vector3r::Zero()
   );

addTranslationRotationConstraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The generic 6-DOF spring-damper joint from which all the specialised joints
above are derived. Each of the six DOF (3 translation, 3 rotation, all in the
joint frame) can be independently stiffened or left free.

.. code-block:: cpp

   size_t idx = engine.addTranslationRotationConstraint(
       entt::entity a,
       entt::entity b,
       const JointFrame& frame,
       Vector3r K_trans = Vector3r::Constant(infinity),  // [Kx, Ky, Kz] N/m
       Vector3r D_trans = Vector3r::Zero(),               // [Dx, Dy, Dz] Ns/m
       Vector3r K_rot   = Vector3r::Zero(),               // [Krx, Kry, Krz] Nm/rad
       Vector3r D_rot   = Vector3r::Zero()                // [Drx, Dry, Drz] Nms/rad
   );

All six stiffness and damping components are expressed in the **joint frame**
(the frame defined by ``frame.A_refJ``):

- ``K_trans[0]`` / ``D_trans[0]`` -- spring/damper along joint x-axis
- ``K_trans[1]`` / ``D_trans[1]`` -- spring/damper along joint y-axis
- ``K_trans[2]`` / ``D_trans[2]`` -- spring/damper along joint z-axis
- ``K_rot[0]``   / ``D_rot[0]``    -- torsion/damper about joint x-axis
- ``K_rot[1]``   / ``D_rot[1]``    -- bending/damper about joint y-axis
- ``K_rot[2]``   / ``D_rot[2]``    -- bending/damper about joint z-axis

Example: a 3-DOF translational spring (compliant joint, rotation free):

.. code-block:: cpp

   engine.addTranslationRotationConstraint(a, b, frame,
       Vector3r(1e3, 1e3, 1e3),  // 1 kN/m in all directions
       Vector3r(10, 10, 10),     // 10 Ns/m damping
       Vector3r::Zero(),         // no rotational stiffness
       Vector3r::Zero()
   );

addBeamConstraint
~~~~~~~~~~~~~~~~~

Connects two beam-segment entities with the Cosserat-rod spring model.
This is called automatically by ``createBeam``; you only need it when
assembling beam chains manually.

.. code-block:: cpp

   size_t idx = engine.addBeamConstraint(
       entt::entity a,
       entt::entity b,
       const BeamSpringParams& springs,
       const BeamCrossSection& section
   );

See :doc:`beams` for full documentation of ``BeamSpringParams`` and
``BeamCrossSection``.

Driving constraint velocities at runtime
-----------------------------------------

Three setters let you prescribe a target velocity for a constraint pattern.
They are typically called inside ``SceneBase::updateScene``.

.. code-block:: cpp

   // Single scalar velocity (for LinearDistanceConstraint)
   engine.setConstraintScalarVelocity(idx, v);

   // Linear velocity in joint frame (for TranslationRotationConstraint)
   engine.setConstraintLinearVelocity(idx, Vector3r(vx, vy, vz));

   // Angular velocity in joint frame
   engine.setConstraintAngularVelocity(idx, Vector3r(wx, wy, wz));

Example: prismatic actuator extending at 0.1 m/s along the joint x-axis:

.. code-block:: cpp

   size_t prismatic = engine.addTranslationalConstraint(base, slider,
       JointFrame::fromAxis(anchorPos, Vector3r::UnitX()));

   // In updateScene:
   engine.setConstraintLinearVelocity(prismatic, Vector3r(0.1, 0, 0));

Quick constraint-type reference
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Method
     - DOF locked
     - Typical use
   * - ``addRigidConstraint``
     - 6 / 6
     - Weld two bodies together.
   * - ``addHingeConstraint``
     - 5 / 6
     - Revolute joint (door hinge, wheel axle, pendulum pin).
   * - ``addTranslationalConstraint``
     - 3 trans / 3
     - Lock position, allow free rotation (ball socket without translation).
   * - ``addRotationConstraint``
     - 3 rot / 3
     - Lock orientation, allow free translation (prismatic guide).
   * - ``addTranslationRotationConstraint``
     - configurable
     - Generic spring-damper joint with full per-axis control.
   * - ``addLinearDistanceConstraint``
     - 1 (distance)
     - Inextensible cable or spring between two attachment points.
   * - ``addBeamConstraint``
     - 6 (compliant)
     - Cosserat-rod link between two beam segments.
