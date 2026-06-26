Trajectories and Kinematic Control
=================================

A trajectory lets you drive an entity from a callback instead of from the contact/constraint solve. The engine stores a :cpp:struct:`C_StaticTrajectory <cardillo::C_StaticTrajectory>` component and updates the entity once per step in the trajectory synchronization code.

.. important::
Calling :cpp:func:`PhysicsEngine::addTrajectory <cardillo::physics::PhysicsEngine::addTrajectory>` on a dynamic body first makes that body static by removing its dynamic mass/inertia tags. The body can still remain collidable and visible, but gravity and solver-driven motion no longer act on it.

Trajectory callback types
-------------------------

The callback aliases are defined in the ECS headers; the public typedefs are
:cpp:typedef:`TrajectoryPose <cardillo::TrajectoryPose>` and
:cpp:typedef:`TrajectoryTwist <cardillo::TrajectoryTwist>`.

.. code-block:: cpp

    using TrajectoryPose  = std::pair<Vector3r, :cpp:type:`Quaternion4r <cardillo::Quaternion4r>`>;
    using TrajectoryTwist = std::pair<Vector3r, Vector3r>;

``TrajectoryPose`` is ``(position, orientation)``.
``TrajectoryTwist`` is ``(linearVelocity, angularVelocity)`` where linear
velocity is in inertial frame and angular velocity follows the engine's usual
Body Basis convention.

How the engine uses them
------------------------

You can provide a pose callback, a twist callback, or both:

- **Pose only** -- the engine samples the pose at ``t`` and ``t + dt`` and
  differentiates it to obtain velocity.
- **Twist only** -- the engine sets the entity's linear and angular velocity,
  then the regular position update advances the pose from that velocity.
- **Pose + twist** -- the engine applies the pose directly and also uses the
  supplied velocity.

If the entity has no orientation or angular-velocity component, the rotational
part is ignored.

Callback trajectories
---------------------

.. code-block:: cpp

   void engine.addTrajectory(
       entt::entity e,
       std::optional<std::function<TrajectoryPose(real_t)>> positionFunc,
       std::optional<std::function<TrajectoryTwist(real_t)>> velocityFunc
   );

Examples:

.. code-block:: cpp

   using namespace cardillo::physics;

   // Position-driven motion
   engine.addTrajectory(
       body,
       std::make_optional<std::function<TrajectoryPose(real_t)>>(
           [](real_t t) -> TrajectoryPose {
               return {
                   Vector3r(0.0, 0.0, 1.0 + 0.1 * std::sin(t)),
                   Quaternion4r::Identity()
               };
           }),
       std::nullopt);

   // Velocity-driven motion
   engine.addTrajectory(
       body,
       std::nullopt,
       std::make_optional<std::function<TrajectoryTwist(real_t)>>(
           [](real_t) -> TrajectoryTwist {
               return {
                   Vector3r(0.5, 0.0, 0.0),
                   Vector3r::Zero()
               };
           }));

   // Pose + velocity
   engine.addTrajectory(
       body,
       std::make_optional<std::function<TrajectoryPose(real_t)>>(
           [](real_t t) -> TrajectoryPose {
               return {
                   Vector3r(std::cos(t), std::sin(t), 0.0),
                   Quaternion4r::Identity()
               };
           }),
       std::make_optional<std::function<TrajectoryTwist(real_t)>>(
           [](real_t t) -> TrajectoryTwist {
               return {
                   Vector3r(-std::sin(t), std::cos(t), 0.0),
                   Vector3r::Zero()
               };
           }));

Spline trajectories
-------------------

There is also a convenience overload for spline-based looping motion:

.. code-block:: cpp

   template <class TSpline>
   void engine.addTrajectory(entt::entity e, const TSpline& spline, real_t period);

This overload:

- accepts spline types derived from ``misc::SplinePattern``
- maps time to ``phase = fmod(max(t, 0), period) / period``
- samples the spline position
- sets orientation to ``Quaternion4r::Identity()``

Example:

.. code-block:: cpp

   cardillo::misc::CatmullRomSpline path;
   path.addControlPoint(Vector3r(0, 0, 0));
   path.addControlPoint(Vector3r(1, 0, 0));
   path.addControlPoint(Vector3r(1, 1, 0));
   path.addControlPoint(Vector3r(0, 1, 0));

   engine.addTrajectory(body, path, 4.0);

Removing a trajectory
---------------------

.. code-block:: cpp

   engine.removeTrajectory(body);

This only removes the :cpp:struct:`C_StaticTrajectory` component. It does **not** recreate
:cpp:struct:`C_PhysicsObject`, mass, or inertia for a body that was made static when the
trajectory was attached.

In practice that means:

- the entity stops following the scripted motion
- its current pose/velocity remain as they are
- it stays static unless you recreate the body or manually restore the missing
  ECS components

Typical uses
------------

Use trajectories for objects that should behave like kinematic scene elements:

- conveyor belt slats
- moving platforms
- scripted grippers or fixtures
- guide geometry for contact tests

If you want a body to remain fully dynamic, do not attach a trajectory to it.
