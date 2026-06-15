Bodies
======

All body-creation methods live on ``PhysicsEngine`` and delegate to the factory
classes in ``src/physics/assets/``. Rigid bodies, static bodies, point masses,
and height fields all return an ``entt::entity`` that you keep around to pass
to constraint factories, trajectory setters, and state queries. 

.. contents:: On this page
   :local:
   :depth: 2

RigidState (initial kinematics)
---------------------------------

``RigidState`` (``src/physics/api/physics_types.hpp``) describes where a body
starts in the world and how fast it is moving:

.. list-table::
   :header-rows: 1
   :widths: 25 50 55

   * - Field
     - Default
     - Meaning
   * - ``position``
     - ``(0, 0, 0)``
     - World-frame centre-of-mass position in metres.
   * - ``orientation``
     - ``Identity``
     - World-frame orientation as a unit quaternion
       (``Eigen::Quaternion<real_t>``).
   * - ``rotation``
     - ``Identity``
     - 3x3 rotation matrix derived from ``orientation``. The named constructors
       keep this in sync automatically; if you assign ``orientation`` directly
       you must update ``rotation`` as well (or prefer the constructors).
   * - ``linearVelocity``
     - ``(0, 0, 0)``
     - World-frame linear velocity in m/s.
   * - ``angularVelocity``
     - ``(0, 0, 0)``
     - Body-frame angular velocity in rad/s.

``RigidState`` provides several convenience constructors so you rarely need
to set every field manually:

.. code-block:: cpp

   using namespace cardillo::physics;

   // Position only
   RigidState s1(Vector3r(0, 0, 1));

   // Position + orientation
   RigidState s2(Vector3r(0, 0, 1),
                 Quaternion4r(AngleAxis3r(M_PI / 4, Vector3r::UnitZ())));

   // Position + linear velocity + orientation
   RigidState s3(Vector3r(0, 0, 1),
                 Vector3r(0, 0, -0.5),
                 Quaternion4r::Identity());

   // Full specification: position, linear velocity, orientation, angular velocity
   RigidState s4(Vector3r(0, 0, 1),
                 Vector3r(0, 0, -0.5),
                 Quaternion4r::Identity(),
                 Vector3r(0, 1, 0));

   // Expressed relative to another entity's frame
   RigidState s5(localPos, localVel, localQuat, localOmega, parentEntity, registry);

.. note::
   For mesh-based bodies the factory recomputes the world position and
   orientation from the mesh's principal-axes frame. The values you pass in
   ``RigidState`` are treated as an offset from the model origin, not the
   centre of mass.

RigidProps (mass and flags)
-----------------------------

``RigidProps`` (``src/physics/api/physics_types.hpp``) controls mass, friction,
and pipeline participation:

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Field
     - Default
     - Meaning
   * - ``mass``
     - ``nullopt``
     - Explicit mass in kg. Takes priority over ``density``. If both
       ``mass`` and ``density`` are unset the body is created with zero mass
       (treated as static by the integrator).
   * - ``density``
     - ``nullopt``
     - Density in kg/m³. The factory multiplies this by the shape's computed
       volume to obtain the mass. Use ``RigidProps::withDensity(rho)`` as a
       named constructor.
   * - ``friction``
     - ``-1``
     - Per-body Coulomb coefficient. A negative value means "use the global
       default" from ``Config::friction_default_mu``. Set a positive value
       to override per body.
   * - ``collidable``
     - ``true``
     - Register this body with the collision system. Set ``false`` for
       invisible constraint anchors or bodies that must never touch anything.
   * - ``visual``
     - ``true``
     - Include this body in VTK output. Set ``false`` to hide helper bodies.

.. code-block:: cpp

   // Explicit mass
   RigidProps p1(1.0);                       // 1 kg, default friction
   RigidProps p2(1.0, 0.3);                  // 1 kg, μ = 0.3
   RigidProps p3(1.0, 0.3, /*vis=*/true, /*coll=*/false); // no collision

   // Density-based
   RigidProps p4 = RigidProps::withDensity(7750.0); // steel

Shape types
-----------

The second argument to ``addRigidBody`` / ``addStaticBody`` is a ``RigidShape``
variant. All shapes are defined in ``src/physics/api/physics_types.hpp``.

CubeShape
~~~~~~~~~

An axis-aligned box. The collision geometry and visual quad both follow the
body's orientation quaternion.

.. code-block:: cpp

   // Half-extents: the box extends ±hx along x, ±hy along y, ±hz along z
   CubeShape(Vector3r halfExtents)

   CubeShape cube(Vector3r(0.5, 0.25, 0.1));  // 1 m × 0.5 m × 0.2 m box

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``halfExtents``
     - Half-extents along local x, y, z in metres.

Inertia is computed as a solid box:
:math:`I_{xx} = m(h_y^2 + h_z^2)/3`, etc.

SphereShape
~~~~~~~~~~~

.. code-block:: cpp

   SphereShape(real_t radius)

   SphereShape s(0.05);  // 5 cm radius sphere

Inertia: :math:`I = \frac{2}{5} m r^2` along all axes.

CapsuleShape
~~~~~~~~~~~~

A cylinder capped with two hemispheres. The long axis runs along the body's
local **z-axis**. Rotate the body's ``orientation`` to align it in world space.

.. code-block:: cpp

   CapsuleShape(real_t radius, real_t halfLength)

   // radius = 0.05 m (cap radius), halfLength = 0.45 m
   // Total length = 2*(halfLength + radius) = 1.0 m
   CapsuleShape arm(0.05, 0.45);

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``radius``
     - Hemisphere cap radius in metres.
   * - ``halfLength``
     - Half-length of the cylindrical shaft between the two caps (metres).

The full extent of the capsule along its local z-axis is
``-(halfLength + radius)`` to ``+(halfLength + radius)``.

Inertia is computed via COAL's ``coal::Capsule::computeMomentofInertia()``,
which accounts for both the cylindrical shaft and the two end caps.

CylinderShape
~~~~~~~~~~~~~

A flat-ended cylinder. Long axis along the body's local **z-axis**.

.. code-block:: cpp

   CylinderShape(real_t radius, real_t halfLength)

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``radius``
     - Barrel radius in metres.
   * - ``halfLength``
     - Half of the total cylinder height in metres.

Inertia:
:math:`I_{zz} = \frac{1}{2}mr^2`,
:math:`I_{xx} = I_{yy} = \frac{m(3r^2 + (2h)^2)}{12}`.

ConeShape
~~~~~~~~~

A right circular cone. Tip points along the body's local **+z** axis.

.. code-block:: cpp

   ConeShape(real_t radius, real_t height)

   ConeShape c(0.1, 0.3);  // base radius 10 cm, height 30 cm

PlaneShape
~~~~~~~~~~

An infinite (collision-wise) flat surface. The visual quad has finite size.
Typically used with ``addStaticBody``.

.. code-block:: cpp

   PlaneShape(Vector3r normal, Vector3r up, real_t sizeX, real_t sizeY)

   // Horizontal floor in the xy-plane, facing +z
   PlaneShape floor(Vector3r::UnitZ(), Vector3r::UnitY(), 5.0, 5.0);

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``normal``
     - World-space outward normal direction (the "above" side of the plane).
   * - ``up``
     - Up direction for the visual quad, used to orient the rendered mesh.
   * - ``sizeX``
     - Visual half-extent along the tangent axis orthogonal to ``normal``
       and ``up``.
   * - ``sizeY``
     - Visual half-extent along ``up``.

.. important::
   ``PlaneShape`` has no mass or inertia. Passing it to ``addRigidBody``
   with a non-zero mass will produce an undefined inertia. Always use
   ``addStaticBody`` for floor planes.

MeshShape
~~~~~~~~~

A triangle mesh loaded from an OBJ or STL file. The engine normalises the
mesh to its principal axes and volume-weighted centre of mass before
computing inertia. Mesh assets are cached so the same path can be used by
multiple bodies without re-loading.

.. code-block:: cpp

   MeshShape(std::string path,
             Vector3r scale = {1,1,1},
             bool use_bbox_collider = false,
             bool show_collider     = false)

   MeshShape link("res/meshes/double_pendulum/link1.stl",
                  Vector3r(1, 1, 1));          // unit scale

   MeshShape approx("res/meshes/complex.obj",
                    Vector3r(0.001, 0.001, 0.001),  // STL in mm → m
                    true,   // use bounding-box collider (faster)
                    false); // don't render the bounding box

.. list-table::
   :header-rows: 1

   * - Field
     - Meaning
   * - ``path``
     - Path to the OBJ or STL file on disk.
   * - ``scale``
     - Per-axis scale applied to the mesh on load. Use ``(0.001, 0.001, 0.001)``
       to convert millimetre STL exports to metres.
   * - ``use_bbox_collider``
     - Replace the exact mesh hull with its axis-aligned bounding box for
       collision. Much faster for complex meshes where exact geometry is not
       needed.
   * - ``show_collider``
     - When using ``use_bbox_collider``, also draw the bounding box in the VTK
       output.

.. note::
   For dynamic mesh bodies, ``RigidState::position`` is interpreted as an
   offset from the model's mesh origin, and the actual world position/orientation
   are adjusted by the mesh's computed centre-of-mass offset and principal-axes
   rotation. Inspect ``debug.mesh = true`` in the config to print these values.

Dynamic vs. static bodies
--------------------------

``addRigidBody``
~~~~~~~~~~~~~~~~

Creates a body that becomes **dynamic** only when the resolved mass is
positive. If mass and density are both unset (or resolve to zero), the factory
still creates the body's pose, shape, visual, and collision components, but it
does not add ``C_PhysicsObject`` / rigid-body mass-inertia data, so the solver
treats the body as static.

.. code-block:: cpp

   entt::entity body = engine.addRigidBody(shape, state, props);

``addStaticBody``
~~~~~~~~~~~~~~~~~

Creates a **static** obstacle with zero mass. It participates in collision
detection but never moves under forces or constraints. Internally this calls
``addRigidBody`` with ``RigidProps(0)``.

.. code-block:: cpp

   entt::entity floor = engine.addStaticBody(
       PlaneShape(Vector3r::UnitZ(), Vector3r::UnitY(), 10, 10),
       RigidState{});  // at origin, default orientation

Other body types
----------------

Point masses
~~~~~~~~~~~~

A point particle with translational (but no rotational) degrees of freedom.
Useful for simple particle simulations or as massless probe points.

.. code-block:: cpp

   entt::entity point = engine.addPointMass(
       real_t mass,        // kg
       Vector3r x0,        // initial position
       Vector3r v0,        // initial velocity
       real_t radius = 0.05  // visual sphere radius
   );

Height-field obstacles
~~~~~~~~~~~~~~~~~~~~~~

A static terrain surface loaded from an OpenEXR file.

.. code-block:: cpp

   entt::entity hf = engine.addObstacleHeightField(
       Vector3r position,         // world position of the field's origin
       Quaternion4r orientation,  // orientation of the field
       std::string exrPath,       // path to the .exr file
       real_t x_dim,              // physical width of the field (metres)
       real_t y_dim,              // physical depth of the field (metres)
       real_t z_scale = 1.0,      // vertical scale factor
       real_t min_height = 0.0    // clip heights below this value
   );

Soft bodies
~~~~~~~~~~~

A deformable body assembled from a surface mesh. Each triangle vertex becomes
a point mass; edges become spring constraints.

.. code-block:: cpp

   std::vector<entt::entity> nodes = engine.addSoftBody(
       std::string objPath,                           // OBJ surface mesh
       real_t stiffness,                              // spring stiffness (N/m)
       real_t damping,                                // spring damping (Ns/m)
       Vector3r position        = Vector3r::Zero(),
       Quaternion4r orientation = Quaternion4r::Identity(),
       Vector3r linearVelocity  = Vector3r::Zero(),
       Vector3r angularVelocity = Vector3r::Zero(),
       real_t totalMass  = 0.0,   // 0 → use per-node defaults
       real_t nodeRadius = 0.02   // visual/collision radius per node
   );

The returned vector contains one ``entt::entity`` per mesh vertex. You can
attach constraints or trajectories to individual nodes.

Disabling collisions between a pair
------------------------------------

After adding bodies, suppress specific pairs that you know will never (or
should never) collide, most commonly a body and its constraint anchor:

.. code-block:: cpp

   engine.disableCollisionBetween(anchor, arm);
