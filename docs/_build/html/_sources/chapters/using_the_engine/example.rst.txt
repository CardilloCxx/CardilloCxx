Worked Example: A Single Pendulum
=================================

This example builds a simple pendulum scene with:

- one dynamic capsule arm
- one small static pivot body
- one hinge joint between pivot and arm
- gravity in ``-z``
- a 30° initial offset
- a small initial angular velocity
- a simple step loop

Using a real pivot body keeps the hinge inside the jointed pair and matches the
common pattern used throughout the example scenes.

Minimal config
--------------

You can load these settings from a config file or assign them in code:

.. code-block:: ini

   sim.T = 5.0
   sim.dt = 1e-3
   sim.gravity = 0 0 -9.81

   output.folder = ./vtk_out
   output.filename_prefix = pendulum
   output.interval_steps = 10

Code
----

.. code-block:: cpp

   #include "config/config.hpp"
   #include "physics/api/physics.hpp"
   #include <Eigen/Geometry>
   #include <cmath>

   using namespace cardillo;
   using namespace cardillo::physics;

   int main() {
       config::Config cfg;
       cfg.sim_T = 5.0;
       cfg.sim_dt = 1e-3;
       cfg.sim_gravity = Vector3r(0, 0, -9.81);
       cfg.output_folder = "./vtk_out";
       cfg.output_filename_prefix = "pendulum";
       cfg.output_interval_steps = 10;

       PhysicsEngine engine(cfg);

       const real_t radius = 0.05;
       const real_t halfLength = 0.45;
       const real_t pivotToCom = radius + halfLength;   // 0.5 m
       const real_t angle = 30.0 * M_PI / 180.0;
       const real_t initialOmega = 0.5;                 // rad/s about hinge axis

       const entt::entity pivot = engine.addStaticBody(
           CubeShape(Vector3r(0.03, 0.03, 0.03)),
           RigidState{});

       const Vector3r position(
           pivotToCom * std::sin(angle),
           0.0,
          -pivotToCom * std::cos(angle));

       const Quaternion4r orientation(
           Eigen::AngleAxis<real_t>(-angle, Vector3r::UnitY()));

       const RigidState armState(
           position,
           Vector3r::Zero(),
           orientation,
           Vector3r(0.0, initialOmega, 0.0));

       const entt::entity arm = engine.addRigidBody(
           CapsuleShape(radius, halfLength),
           armState,
           RigidProps(1.0));

       engine.addHingeConstraint(
           pivot,
           arm,
           JointFrame::fromAxis(Vector3r::Zero(), Vector3r::UnitY()));

       engine.disableCollisionBetween(pivot, arm);
       engine.track(arm, "pendulum_arm");

       while (!engine.isFinished()) {
           engine.step();
       }

       return 0;
   }


See also: existing scenes
-------------------------

For more complete setups, take a look at the scenes already in
``examples/scenes``.

Those examples show real project usage of hinges, trajectories, beams, and
larger assembled mechanisms.
