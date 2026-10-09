Trajectories and Kinematic Control
=================================

A trajectory lets you drive an entity from a callback instead of from the contact/constraint solve. The engine stores a :cpp:struct:`C_StaticTrajectory <cardillo::C_StaticTrajectory>` component and updates the entity once per step in the trajectory synchronization code.

.. important::
Calling :cpp:func:`PhysicsEngine::addTrajectory <cardillo::physics::PhysicsEngine::addTrajectory>` on a dynamic body first makes that body static by removing its dynamic mass/inertia tags. The body can still remain collidable and visible, but gravity and solver-driven motion no longer act on it.

Trajectory callback type
------------------------

A trajectory prescribes the pose of the entity as a function of time,
:cpp:typedef:`TrajectoryPose <cardillo::TrajectoryPose>`:

.. code-block:: cpp

    using TrajectoryPose = std::pair<Vector3r, Quaternion4r>;  // (position, orientation)

How the engine uses it
----------------------

Only the pose is prescribed; the velocity is always derived from it, so that
pose and velocity are consistent. In every Moreau step
:math:`t_n \to t_n + \Delta t`, the engine

- sets the pose to its value at the intermediate time
  :math:`t_{n+\theta} = t_n + (1-\theta)\,\Delta t`, where all position-level
  measures of the step are evaluated, and
- sets the velocity to the mean velocity over the step,
  ``(pose(t_n + dt) - pose(t_n)) / dt`` (angular velocity in the body frame,
  from the relative rotation of the two orientations).

The velocity enters the force laws of constraints and contacts with the driven
body as a known source term. Velocity-only trajectories are not supported:
prescribe a velocity by its integrated pose, e.g. a velocity
``v(t) = v0 sin(w t)`` by the position ``x0 + v0 / w (1 - cos(w t))``, and a
constant angular velocity ``w`` about the body x-axis by the orientation
``q0 * AngleAxis(w t, UnitX)``. If the entity has no orientation or
angular-velocity component, the rotational part is ignored.

Callback trajectories
---------------------

.. code-block:: cpp

   void engine.addTrajectory(entt::entity e, std::function<TrajectoryPose(real_t)> positionFunc);

Examples:

.. code-block:: cpp

   using namespace cardillo::physics;

   // Vertical oscillation
   engine.addTrajectory(body, [](real_t t) -> TrajectoryPose {
       return {Vector3r(0.0, 0.0, 1.0 + 0.1 * std::sin(t)), Quaternion4r::Identity()};
   });

   // Circular path with a constant spin of 2 rad/s about the body z-axis
   engine.addTrajectory(body, [](real_t t) -> TrajectoryPose {
       return {Vector3r(std::cos(t), std::sin(t), 0.0),
               Quaternion4r(Eigen::AngleAxis<real_t>(2.0 * t, Vector3r::UnitZ()))};
   });

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
