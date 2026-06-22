Getting Started
===============

Cardillo builds as a C++17 library with example scenes. The fastest path is to configure the root CMake project, build it, and run one of the scene configs under ``examples/scenes``.

Install and build
-----------------

Required on Linux:

- CMake 3.16 or newer
- A C++17 compiler such as GCC or Clang
- Python 3

From the repository root:

.. code-block:: bash

	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	cmake --build build -j$(nproc)

If you want the Clarabel backend, install Rust and Cargo before configuring the build.

.. code-block:: bash

	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Run an example
--------------

Each example scene is driven by a config file that sets the scene name, output folder, solver, and time step. For example:

.. code-block:: bash

	./build/bin/main ./examples/scenes/wilberforce/scene.config

If you want a compact code walkthrough after building the project, see the
worked pendulum example in :doc:`using_the_engine/example`.

Docs build
----------

The documentation is generated with Sphinx, Breathe, and Exhale. Rebuild it from ``docs/`` after changing the public headers or chapter sources:

.. code-block:: bash

	uv run make html

For details on all types (vectors, matrices, shapes, components), see :doc:`types`.
