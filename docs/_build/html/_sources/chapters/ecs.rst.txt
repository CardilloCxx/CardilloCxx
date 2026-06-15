Entity Component System
=======================

Cardillo uses EnTT to describe every simulation object as a plain entity
(``entt::entity``) carrying only the data it needs. Instead of classes or inheritance, just components and tags are queried by systems. All public component definitions live in
``src/physics/ecs_types.hpp``.

.. contents:: On this page
   :local:
   :depth: 2

The ECS mental model
---------------------

A system does not own objects; it **reads** the combination of components that
identify a behaviour:

- ``C_RigidBodyTag + C_Mass + C_PhysicsObject`` → rigid body, participate in dynamics.
- ``C_Collidable`` → include in collision detection.
- ``C_VisualObject`` → export to VTK output.
- ``C_BeamElement`` → connect into a beam chain.

The engine's own systems (collision broad-phase, solver assembly, integrator,
VTK writer) are all just views over the registry.

Composing an entity: rigid body factory snippet
------------------------------------------------

Here is how ``RigidBodyFactory::create`` (``src/physics/assets/rigid_body_factory.cpp``)
builds a dynamic capsule from scratch:

.. code-block:: cpp

   // RigidBodyFactory::create in rigid_body_factory.cpp
   auto& reg = system.ecs();
   const auto e = reg.create();                       // 1. bare entity

   // 2. Core kinematic state (every rigid body gets these)
   reg.emplace<C_Position3>(e, state.position);       // world pos (m)
   reg.emplace<C_Orientation>(e, state.orientation);  // unit quaternion
   reg.emplace<C_LinearVelocity3>(e, state.linearVelocity);
   reg.emplace<C_AngularVelocity3>(e, state.angularVelocity);

   // 3. Pipeline participation tags
   if (props.visual)        reg.emplace<C_VisualObject>(e);     // VTK output
   if (props.collidable)    reg.emplace<C_Collidable>(e);       // collision detector

   const real_t mu = (props.friction >= 0) ? props.friction : cfg.friction_default_mu;
   reg.emplace<C_Friction>(e, mu);                           // Coulomb coefficient

   // 4. Dynamic body components (only when mass > 0)
   reg.emplace<C_PhysicsObject>(e);
   reg.emplace<C_RigidBodyTag>(e);
   reg.emplace<C_Mass>(e, mass);
   reg.emplace<C_InertiaDiag>(e, inertiaDiag);

   // 5. Shape-specific components (inside the std::visit<Visitor> dispatch)
   // For a CapsuleShape:
   reg.emplace<C_Capsule>(e, shape.radius, shape.halfLength);
   reg.emplace<C_RB_Capsule>(e, shape.radius, shape.halfLength);
   reg.emplace<C_CapsuleVisualTag>(e);

The key takeaway: a rigid body entity carries roughly 14–15 components. A system
that needs to read only pose and friction would query exactly those two:

.. code-block:: cpp

   engine.ecs().view<const C_Position3, const C_Friction>().each([](auto e,
         const C_Position3& p, const C_Friction& f) { /* ... */ });

All component types
-------------------

Below every public component is listed with its data members. Tags (empty structs)
are grouped together; data components show their fields and units.

Data components
~~~~~~~~~~~~~~~

Position and orientation
""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_Position3``
     - ``Vector3r value``
     - World-space position (m).
   * - ``C_Orientation``
     - ``Quaternion4r value``
     - World-space orientation as a unit quaternion.
   * - ``C_DirectorTriad``
     - ``Matrix33r value``
     - Body-frame triad (used for director-type bodies).

Velocity and acceleration
"""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_LinearVelocity3``
     - ``Vector3r value``
     - World-space linear velocity (m/s).
   * - ``C_AngularVelocity3``
     - ``Vector3r value``
     - Body-frame angular velocity (rad/s).
   * - ``C_LinearAcceleration3``
     - ``Vector3r value``
     - Current linear acceleration accumulator.
   * - ``C_AngularAcceleration3``
     - ``Vector3r value``
     - Current angular acceleration accumulator.

Mass and inertia
"""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_Mass``
     - ``real_t m``
     - Body mass (kg). Zero → static body.
   * - ``C_InertiaDiag``
     - ``Vector3r I``
     - Body-frame diagonal inertia: :math:`[I_{xx},\, I_{yy},\, I_{zz}]`.

Shape data components
"""""""""""""""""""""

Rigid-body primitives (local geometry in the body frame):

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Component
     - Fields
     - Meaning
   * - ``C_Cube``
     - ``Vector3r center``, ``halfExtents``, ``Quaternion4r q``
     - Box geometry + pose.
   * - ``C_Capsule``
     - ``real_t radius``, ``halfLength``
     - Capsule shape (z-axis long).
   * - ``C_Cylinder``
     - ``real_t radius``, ``halfLength``
     - Cylinder shape (z-axis long).
   * - ``C_Cone``
     - ``real_t radius``, ``height``
     - Cone shape. Tip along +z.
   * - ``C_Plane``
     - ``Vector3r normal``, ``up``, ``sizeX``, ``sizeY``
     - Plane surface for visualisation.

Rigid-body collision primitives (used by the collision manager):

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Component
     - Fields
     - Meaning
   * - ``C_RB_Cube``
     - ``Vector3r center``, ``halfExtents``, ``q``
     - AABB-aligned cube collider.
   * - ``C_RB_Plane``
     - ``normal``, ``up``, ``sizeX``, ``sizeY``
     - Plane collider.
   * - ``C_RB_Capsule``
     - ``real_t radius``, ``halfLength``
     - Capsule collider (via COAL).
   * - ``C_RB_Cylinder``
     - ``real_t radius``, ``halfLength``
     - Cylinder collider via COAL.
   * - ``C_RB_Cone``
     - ``real_t radius``, ``height``
     - Cone collider via COAL.
   * - ``C_RB_Sphere``
     - *(none)*
     - Marker: this body has a sphere collider.
   * - ``C_RB_Mesh``
     - *(none)*
     - Marker: this body has a mesh (BVH) collider.

Other shape-related data components:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Component
     - Fields
     - Meaning
   * - ``C_Radius``
     - ``real_t r``
     - Visual/collision sphere radius (used by capsule visualisation).
   * - ``C_Mesh``
     - ``std::string path``, ``Vector3r scale``
     - Mesh asset path and per-axis scale.
   * - ``C_HeightField``
     - ``path``, ``x_dim``, ``y_dim``, ``z_scale``, ``min_height``
     - Height-field terrain parameters.
   * - ``C_SoftBodySurface``
     - ``triangles`` (``Vector3i[]``), ``nodes`` (``entity[]``)
     - Surface mesh for a soft body: one triangle per face, one node per vertex.

Physical property data components
"""""""""""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_Friction``
     - ``real_t mu``
     - Coulomb friction coefficient for this body.

Trajectory data component
"""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_StaticTrajectory``
     - ``positionFunc``, ``velocityFunc``, ``elapsed``, ``initialized``, ``previousPosition``
     - Kinematic override. When present the integrator writes pose/velocity from the function instead of solving forces.

Force/torque data components
""""""""""""""""""""""""""""

These are written by external callers and read by the Moreau integrator each step:

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_ExternalForce``
     - ``Vector3r f``
     - World-frame force applied this step (N).
   * - ``C_ExternalTorque``
     - ``Vector3r tau``
     - World-frame torque applied this step (Nm).

Output and bookkeeping data components
"""""""""""""""""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_TrackTag``
     - ``std::string name``
     - Label for CSV track output. Set via ``engine.track(entity, name)``.
   * - ``C_BodyIndex``
     - ``int b``
     - Body index in the assembled mass matrix (set during pipeline assembly).

Beam component
"""""""""""""""""""""""""""""""""""

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Component
     - Fields
     - Meaning
   * - ``C_BeamElement``
     - ``prev``, ``next`` (optional entities), ``l0`` (rest length), ``l`` (current length)
     - Links this body into a beam chain. Filled by the beam factory.

Tags (presence indicates behaviour)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Tag
     - Meaning
   * - ``C_PhysicsObject``
     - This entity participates in the physics pipeline (dynamics solver). Zero-mass bodies do not get this tag.
   * - ``C_RigidBodyTag``
     - This entity is a rigid body (6 DOF). Together with ``C_Mass`` it identifies dynamic bodies.
   * - ``C_PointMassTag``
     - This entity is a point particle (3 translational DOF only, no rotation).
   * - ``C_RigidBodyDirectorTag``
     - Director-type rigid body (supports additional orientation degrees of freedom).
   * - ``C_Collidable``
     - Include this body in collision detection. Removing it makes the body invisible to the narrow-phase without touching other components.
   * - ``C_VisualObject``
     - Render / export this body in VTK output.
   * - ``C_HeightFieldVisualTag``
     - Mark a height-field for visual rendering.
   * - ``C_RB_HeightField``
     - Mark this entity as a rigid-body height-field (immovable terrain).
   * - ``C_SoftBodyVisualTag``
     - Soft-body node should be rendered as a sphere.

Shape-specific visual tags
"""""""""""""""""""""""""""

Each primitive shape has a dedicated visual tag so systems can branch on body type:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Tag
     - Meaning
   * - ``C_PointVisualTag``
     - Sphere visual (used for sphere shapes).
   * - ``C_CubeVisualTag``
     - Box visual quad.
   * - ``C_CapsuleVisualTag``
     - Capsule mesh visual.
   * - ``C_CylinderVisualTag``
     - Cylinder mesh visual.
   * - ``C_ConeVisualTag``
     - Cone mesh visual.
   * - ``C_PlaneVisualTag``
     - Plane visual quad.
   * - ``C_MeshVisualTag``
     - Loaded triangle-mesh visual.

Querying and iterating
----------------------

The most common EnTT patterns in Cardillo are **views** (filtered reads) and
**storage writes**.

Filtered query to iterate all dynamic rigid bodies:

.. code-block:: cpp

   auto view = engine.ecs().view<
       const C_Mass,
       const C_Position3,
       const C_Orientation
   >();
   // Only entities carrying ALL three components
   for (auto entity : view) {
       real_t mass  = view.get<C_Mass>(entity).m;
       auto   pos   = view.get<C_Position3>(entity).value.transpose();
       auto   quat  = view.get<C_Orientation>(entity).value;
   }

Filtered query to iterate only collidable bodies:

.. code-block:: cpp

   auto collision_view = engine.ecs().view<const C_Collidable, const C_Mass>();

Query with "any of": find all beam segments regardless of their body type:

.. code-block:: cpp

   engine.ecs().any_of<C_BeamElement>(entity);    // boolean check

Adding and removing components at runtime
-----------------------------------------

You can extend an entity's behaviour mid-simulation by emplacing or removing
components directly through the registry. The engine's dirty flags keep assembly
in sync:

.. code-block:: cpp

   auto& ecs = engine.ecs();

   // Give a floor body dynamic behaviour
   ecs.get_or_emplace<C_Mass>(floorEntity, 10.0);
   ecs.emplace<C_PhysicsObject>(floorEntity);
   ecs.emplace<C_RigidBodyTag>(floorEntity);
   engine.markStructureDirty();            // next step will rebuild

   // Remove collision from an entity
   ecs.remove<C_Collidable>(entity);

.. caution::
   Manually removing ``C_Mass`` (setting it to 0) or removing
   ``C_PhysicsObject`` + ``C_RigidBodyTag`` makes the integrator treat the
   entity as a static obstacle.

Scene composition reference
----------------------------

Every scene in ``examples/scenes/`` is an entity composition recipe. As a concrete
example, here is what the **pendulum** scene from :doc:`using_the_engine/example` produces:

.. code-block:: text

   Entity 1 (pivot, static cube):
     C_Cube,        ← visual box geometry
     C_RB_Cube,     ← collider box geometry
     C_VisualObject ← write to VTK
     C_Collidable  ← participate in collision

   Entity 2 (arm, dynamic capsule):
     C_Position3         (0.25, 0, -0.433)
     C_Orientation       (q from -30° around Y)
     C_LinearVelocity3   (0, 0, 0)
     C_AngularVelocity3  (0, 0, 0)
     C_LinearAcceleration3 (zero init)
     C_AngularAcceleration3 (zero init)
     C_Mass              (1.0 kg)
     C_InertiaDiag       (computed from capsule volume)
     C_PhysicsObject     ← dynamics solver
     C_RigidBodyTag      ← rigid body
     C_Capsule           ← shape tag
     C_RB_Capsule        ← collider tag
     C_CapsuleVisualTag  ← render as capsule mesh
     C_Collidable        ← collision detection
     C_Friction          (mu from config)
     C_TrackTag          ("pendulum_arm")

The two entities are linked by a ``TranslationRotationConstraint`` pattern with a
hinge axis along Y at the origin, no ECS component needed for constraints, they
live in ``World::m_constraints_new`` as ``unique_ptr<ConstraintPattern>``.
