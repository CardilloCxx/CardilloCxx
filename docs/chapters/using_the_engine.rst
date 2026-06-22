Using The Engine
================

:cpp:class:`PhysicsEngine <cardillo::physics::PhysicsEngine>` is the single public entry point for
building and advancing simulations. It owns the ECS world, the collision manager,
and the physics pipeline, and exposes them through a clear set of factory methods
and runtime controls.

 The mental model is simple: you create entities through factory methods, connect
 them with constraints or trajectories, then advance the whole ECS-based world in
 a fixed-step loop with :cpp:func:`PhysicsEngine::step <cardillo::physics::PhysicsEngine::step>`. The engine handles collision detection,
assembly, solving, and output, while you focus on scene construction and runtime
control.

The chapters below go through each part of the API in detail. Parameter meanings
are derived from the public API headers and documentation. The final chapter
solves a concrete problem (a swinging pendulum) from scratch.

.. toctree::
   :maxdepth: 2

   using_the_engine/setup
   using_the_engine/bodies
   using_the_engine/constraints
   using_the_engine/splines
   using_the_engine/trajectories
   using_the_engine/beams
   using_the_engine/stepping
   using_the_engine/example

For a complete list of all types, shapes, and components, see :doc:`types`.
