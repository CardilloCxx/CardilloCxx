Splines
=======

Splines are reusable geometric paths used by the engine in two main places:

- as input geometry for ``createBeam`` / ``createBeams``
- as looping kinematic paths for ``addTrajectory(entity, spline, period)``

All spline types live in ``src/misc/spline.hpp`` under ``cardillo::misc``.

.. contents:: On this page
   :local:
   :depth: 2

The common interface
--------------------

All spline classes derive from ``misc::SplinePattern``:

.. code-block:: cpp

   class SplinePattern {
   public:
       virtual ~SplinePattern() = default;
       virtual real_t totalLength() const = 0;
       virtual bool isLoop() const = 0;
       virtual SplineSample sample(real_t alpha) const = 0;
       virtual Vector3r centerOfMass() const = 0;
   };

``sample(alpha)`` returns a ``SplineSample``:

.. code-block:: cpp

   struct SplineSample {
       Vector3r position;
       Vector3r tangent;
       Vector3r normal;
       Vector3r binormal;
   };

The parameter ``alpha`` is normalized to ``[0, 1]`` and clamped by the built-in
implementations.

Built-in spline types
---------------------

LinearSpline
~~~~~~~~~~~~

A straight segment between two points.

.. code-block:: cpp

   misc::LinearSpline line(p0, p1);

Constructor:

.. code-block:: cpp

   LinearSpline(const Vector3r& p0, const Vector3r& p1);

Properties:

- ``totalLength()`` is ``|p1 - p0|``
- ``isLoop()`` is always ``false``
- ``centerOfMass()`` is the midpoint

CircleSpline
~~~~~~~~~~~~

A circular arc in a plane, or a full loop if the angle span is ``2π``.

.. code-block:: cpp

   misc::CircleSpline arc(
       Vector3r::Zero(),
       0.5,
       Vector3r::UnitZ(),
       Vector3r::UnitX(),
       0.0,
       M_PI);

Constructor:

.. code-block:: cpp

   CircleSpline(
       const Vector3r& center,
       real_t radius,
       const Vector3r& normal = Vector3r::UnitZ(),
       const Vector3r& dir0 = Vector3r::UnitX(),
       real_t thetaStart = 0,
       real_t thetaSpan = 2 * M_PI);

Notes:

- ``normal`` defines the plane
- ``dir0`` defines the radial direction at ``alpha = 0``
- ``isLoop()`` is true only when ``|thetaSpan|`` is effectively ``2π``
- ``centerOfMass()`` returns the circle center

HelixSpline
~~~~~~~~~~~

A helical path around an axis.

.. code-block:: cpp

   misc::HelixSpline helix(
       Vector3r::Zero(),
       Vector3r::UnitZ(),
       0.02,
       0.01,
       8.0,
       Vector3r::UnitX());

Constructor:

.. code-block:: cpp

   HelixSpline(
       const Vector3r& center,
       const Vector3r& axisDir,
       real_t radius,
       real_t pitch,
       real_t turns,
       const Vector3r& dir0 = Vector3r::UnitX());

Notes:

- ``axisDir`` is normalized internally
- ``pitch`` is the axial advance per revolution
- ``turns`` is the total number of revolutions
- ``isLoop()`` is always ``false``

CatmullRomSpline
~~~~~~~~~~~~~~~~

A smooth curve through a list of control points.

.. code-block:: cpp

   std::vector<Vector3r> pts = {
       Vector3r(0, 0, 0),
       Vector3r(0.3, 0, 0.1),
       Vector3r(0.6, 0.2, 0.1),
       Vector3r(1.0, 0, 0)
   };

   misc::CatmullRomSpline curve(pts, false);

Constructor:

.. code-block:: cpp

   CatmullRomSpline(std::vector<Vector3r> controlPoints, bool loop);

Notes:

- open and closed curves are supported
- the implementation builds an arc-length lookup table, so ``alpha`` is sampled
  more uniformly by length than raw segment parameter space
- ``centerOfMass()`` is approximated numerically from sampled arc segments

Loading splines from BCC
------------------------

The utility below loads one or more Catmull-Rom splines from a ``.bcc`` file:

.. code-block:: cpp

   std::vector<std::shared_ptr<misc::SplinePattern>> splines =
       misc::loadSplinesFromBCC("res/paths/my_curve.bcc", 0.001);

Signature:

.. code-block:: cpp

   std::vector<std::shared_ptr<SplinePattern>>
   loadSplinesFromBCC(const std::string& filePath, real_t scale = 1);

The ``scale`` factor is applied to all loaded control points.

Using splines to create beams
-----------------------------

A single spline can be sampled into a beam chain with ``createBeam``:

.. code-block:: cpp

   misc::LinearSpline spline(Vector3r(0, 0, 0), Vector3r(1, 0, 0));

   auto [root, tip] = engine.createBeam(
       spline,
       section,
       springs,
       RigidState{},
       RigidProps::withDensity(1200.0),
       32);

For multiple connected spline pieces, use ``createBeams``:

.. code-block:: cpp

   misc::LinearSpline part1(Vector3r(0, 0, 0), Vector3r(0.5, 0, 0));
   misc::CircleSpline part2(Vector3r(0.5, 0.1, 0), 0.1,
                            Vector3r::UnitZ(), Vector3r::UnitX(),
                            -M_PI_2, M_PI_2);

   std::vector<const misc::SplinePattern*> parts{&part1, &part2};

   auto [root, tip] = engine.createBeams(
       parts,
       section,
       springs,
       RigidState{},
       props,
       64);

How it works:

- each spline is sampled immediately during the call
- the total segment budget is distributed across the input splines by length
- adjacent generated beam segments are connected with beam constraints
- adjacent beam-segment collisions are disabled automatically

.. note::
   ``createBeams`` takes raw ``SplinePattern`` pointers, but it samples them
   synchronously inside the call. The spline objects only need to stay alive for
   the duration of that call.

Using splines for trajectories
------------------------------

You can also drive an entity along a spline with the spline overload of
``addTrajectory``:

.. code-block:: cpp

   misc::CircleSpline path(Vector3r::Zero(), 0.5);
   engine.addTrajectory(body, path, 4.0);

Signature:

.. code-block:: cpp

   template <class TSpline>
   void engine.addTrajectory(entt::entity e, const TSpline& spline, real_t period);

Behavior of this overload:

- it copies the spline internally
- it maps time to ``phase = fmod(max(t, 0), period) / period``
- it uses ``sample(phase).position`` as the commanded position
- it always sets orientation to ``Quaternion4r::Identity()``
- it is equivalent to a position-only trajectory callback

If you need non-identity orientation along the path, use the callback-based
trajectory API and compute both pose and velocity yourself.
