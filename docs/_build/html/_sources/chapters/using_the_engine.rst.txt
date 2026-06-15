Using The Engine
================

``PhysicsEngine`` (``src/physics/api/physics_engine.hpp``) is the single public entry point for
building and advancing simulations. It owns the ECS world, the collision manager,
and the physics pipeline, and exposes them through a clear set of factory methods
and runtime controls.

The mental model is simple: you create entities through factory methods, connect
them with constraints or trajectories, then advance the whole ECS-based world in
a fixed-step loop with ``step()``. The engine handles collision detection,
assembly, solving, and output, while you focus on scene construction and runtime
control.

The chapters below go through each part of the API in detail. Parameter meanings
are derived directly from the source headers in ``src/physics/api/`` and
``src/physics/constraints/``. The final chapter solves a concrete problem
(a swinging pendulum) from scratch.

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
