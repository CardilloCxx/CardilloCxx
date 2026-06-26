Types Reference
===============

Cardillo builds on a thin layer of type aliases over **Eigen** and the C++ STL. This chapter explains every publicly relevant type and shows where it lives in the codebase. All names link to the generated API reference.

The types are defined in `misc/types.hpp` and
`physics/api/physics_types.hpp`. 

Core scalar and index types
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``real_t``
     - ``double``. The engine's floating-point type everywhere (positions, velocities, forces, configuration values).
   * - ``index_t``
     - ``int``. Used for entity indices, array indices, and loop counters.

Eigen vector aliases
--------------------

Fixed-size vectors (most commonly used):

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Type
     - Dimensions
     - Typical use
   * - :cpp:type:`Vector2r <cardillo::Vector2r>`
     - 2x1
     - Rarely used directly; mostly for image/spline intermediates.
   * - :cpp:type:`Vector3r <cardillo::Vector3r>`
     - 3x1
     - The workhorse: positions, velocities, normals, forces, torques.
   * - :cpp:type:`Vector4r <cardillo::Vector4r>`
     - 4x1
     - Extended state vectors.
   * - :cpp:type:`VectorXr <cardillo::VectorXr>`
     - Dynamic
     - Stack of all body states, force vectors, impulse vectors.

Integer variants:

.. list-table::
   :header-rows: 1
   :widths: 25 45

   * - ``Vector3i``, ``VectorXi``
     - Integer versions for indices and mesh vertex IDs.
   * - ``Array3b``, ``VectorXb``
     - Boolean masks (e.g. collision flags).

Eigen matrix aliases
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 45

   * - Type
     - Dimensions
     - Typical use
   * - :cpp:type:`Matrix33r <cardillo::Matrix33r>`
     - 3x3
     - Inertia tensors, rotation matrices, Jacobian blocks.
   * - :cpp:type:`Matrix44r <cardillo::Matrix44r>`
     - 4x4
     - Homogeneous transforms.
   * - :cpp:type:`MatrixXXr <cardillo::MatrixXXr>`
     - DynamicxDynamic
     - Full mass matrix, assembled constraint matrices.

Skew-symmetric helper:

A 3x3 matrix representing the cross-product operator, constructed via the helper function :cpp:func:`skew_from_vector <cardillo::skew_from_vector>`.

Quaternion and rotation
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`Quaternion4r <cardillo::Quaternion4r>`
     - Unit quaternion for orientation storage and interpolation.
   * - :cpp:type:`AngleAxis3r <cardillo::AngleAxis3r>`
     - Axis-angle representation; convenient for building rotations from an axis and an angle.

Dynamic-size arrays
-------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`ArrayXr <cardillo::ArrayXr>`
     - Dynamic Eigen array (element-wise operations).
   * - :cpp:type:`ArrayXXr <cardillo::ArrayXXr>`
     - Dynamic Eigen matrix for array-style math.
   * - :cpp:type:`VectorXi <cardillo::VectorXi>`
     - Dynamic integer vector.
   * - :cpp:type:`ArrayXi <cardillo::ArrayXi>`
     - Dynamic integer array.

References (pass-by-reference helpers)
--------------------------------------

These are thin wrappers over ``Eigen::Ref`` that let you pass dynamically-sized data to functions without copying:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`RefVectorXr <cardillo::RefVectorXr>`
     - Read-only reference to a dynamic ``VectorXr``.
   * - :cpp:type:`RefVectorXi <cardillo::RefVectorXi>`
     - Read-only reference to a dynamic ``VectorXi``.
   * - :cpp:type:`RefMatrixXXr <cardillo::RefMatrixXXr>`
     - Read-only reference to a dynamic ``MatrixXXr``.

Maps (non-owning views into external memory)
---------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`MapVectorXr <cardillo::MapVectorXr>`
     - Non-owning view onto an externally allocated ``VectorXr`` buffer.
   * - :cpp:type:`MapMatrixXXr <cardillo::MapMatrixXXr>`
     - Non-owning view onto an externally allocated ``MatrixXXr`` buffer.

Sparse matrices
---------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`CscMatrix <cardillo::CscMatrix>`
     - Compressed Sparse Column format (default for the solver).
   * - :cpp:type:`CsrMatrix <cardillo::CsrMatrix>`
     - Compressed Sparse Row format.
   * - :cpp:type:`SparseMatrix <cardillo::SparseMatrix>`
     - Generic sparse matrix template alias.

STL container aliases
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:type:`StlVectorXr <cardillo::StlVectorXr>`
     - ``std::vector<real_t>`` — generic real arrays.
   * - :cpp:type:`StlVectorXi <cardillo::StlVectorXi>`
     - ``std::vector<index_t>`` — generic integer arrays.
   * - :cpp:type:`StlVectorVectorXr <cardillo::StlVectorVectorXr>`
     - ``std::vector<VectorXr>`` — list of vectors (e.g. per-body data).

Physics-specific types (``cardillo::physics``)
-----------------------------------------------

These live in ``physics/api/physics_types.hpp`` and are the primary types users interact with:

Initial kinematics
~~~~~~~~~~~~~~~~~~

The primary initial state type for a rigid body:

* :cpp:struct:`RigidState <cardillo::physics::RigidState>` — describes where a body starts (position, orientation, linear velocity, angular velocity).

Shape types
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:struct:`CubeShape <cardillo::physics::CubeShape>`
     - Axis-aligned box (half-extents).
   * - :cpp:struct:`SphereShape <cardillo::physics::SphereShape>`
     - Sphere.
   * - :cpp:struct:`CapsuleShape <cardillo::physics::CapsuleShape>`
     - Cylinder capped with hemispheres.
   * - :cpp:struct:`CylinderShape <cardillo::physics::CylinderShape>`
     - Flat-ended cylinder.
   * - :cpp:struct:`ConeShape <cardillo::physics::ConeShape>`
     - Right circular cone.
   * - :cpp:struct:`PlaneShape <cardillo::physics::PlaneShape>`
     - Infinite flat surface (visual quad is finite).
   * - :cpp:struct:`MeshShape <cardillo::physics::MeshShape>`
     - Triangle mesh loaded from OBJ/STL.

The union type:

:cpp:type:`RigidShape <cardillo::physics::RigidShape>` is ``std::variant<CubeShape, PlaneShape, CapsuleShape, CylinderShape, ConeShape, SphereShape, MeshShape>`` — a shape can be any of the above.

Body properties
~~~~~~~~~~~~~~~

The primary properties type for a rigid body:

* :cpp:struct:`RigidProps <cardillo::physics::RigidProps>` — mass, density, friction coefficient, collidable and visual flags.

Beam types
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:enum:`BeamBodyType <cardillo::physics::BeamBodyType>`
     - Enum: ``Cube``, ``Capsule``, ``Cylinder`` — controls the collision/visual shape of beam segments.
   * - :cpp:struct:`BeamCrossSection <cardillo::physics::BeamCrossSection>`
     - Cross-section geometry (width, height, type) with derived properties (area, moments of inertia).
   * - :cpp:struct:`BeamSpringParams <cardillo::physics::BeamSpringParams>`
     - Elastic and damping parameters for the Cosserat-rod beam model.

Other types
-----------

Configuration:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:struct:`Config <cardillo::config::Config>`
     - All global simulation parameters (time step, gravity, solver settings, output).
   * - :cpp:class:`ConfigReader <cardillo::config::ConfigReader>`
     - Loads ``Config`` from a text file.

Collision:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:struct:`Contact <cardillo::collision::Contact>`
     - Information for a single pairwise contact (points, normals, impulses).
   * - :cpp:type:`ContactList <cardillo::collision::ContactList>`
     - ``std::vector<Contact>`` — list of all contacts for a step.
   * - :cpp:struct:`SphereCollider <cardillo::collision::SphereCollider>`
     - Sphere collider resolved from ECS components.
   * - :cpp:struct:`PlaneCollider <cardillo::collision::PlaneCollider>`
     - Plane collider resolved from ECS components.
   * - :cpp:struct:`AabbCollider <cardillo::collision::AabbCollider>`
     - AABB collider resolved from ECS components.
   * - :cpp:struct:`ObbCollider <cardillo::collision::ObbCollider>`
     - OBB collider resolved from ECS components.

Splines:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - :cpp:class:`SplinePattern <cardillo::misc::SplinePattern>`
     - Abstract base for all spline types.
   * - :cpp:struct:`SplineSample <cardillo::misc::SplineSample>`
     - Result of sampling a spline: position, tangent, normal, binormal.
   * - :cpp:class:`LinearSpline <cardillo::misc::LinearSpline>`
     - Straight-line spline between two points.
   * - :cpp:class:`CircleSpline <cardillo::misc::CircleSpline>`
     - Circular arc spline.
   * - :cpp:class:`HelixSpline <cardillo::misc::HelixSpline>`
     - Helical path spline.
   * - :cpp:class:`CatmullRomSpline <cardillo::misc::CatmullRomSpline>`
     - Catmull-Rom spline through control points.

ECS components: all public component definitions live in ``ecs_types.hpp``. Each ``C_*`` type is a data component or tag in the ``cardillo`` namespace. For example:

* Data components: :cpp:type:`C_Position3 <cardillo::C_Position3>`, :cpp:type:`C_Mass <cardillo::C_Mass>`, etc.
* Tags: :cpp:type:`C_PhysicsObject <cardillo::C_PhysicsObject>`, :cpp:type:`C_Collidable <cardillo::C_Collidable>`, etc.
