Getting Started
===============

CardilloCxx builds as a C++20 shared library ``libcardillo`` with example scenes and tests. Its
headers are included as ``<cardillo/...>``.

Requirements
------------

- CMake 3.24 or newer and a C++20 compiler (GCC or Clang)
- VTK (a system VTK is used if found, otherwise it is built from source; a Debian/Ubuntu ParaView
  package that ships its own VTK also works)
- Boost, Assimp and Octomap development packages (needed by Coal)
- OpenMP (optional)
- Rust/cargo, only for the optional Clarabel backend

Eigen 5 and EnTT are taken from an installation if CMake finds them, otherwise they are fetched.
Coal is fetched from the fork ``41pha1/coal``, which adds the custom shapes used by the beam hull
colliders. On Ubuntu:

.. code-block:: bash

	sudo apt install build-essential cmake git libvtk9-dev libboost-dev libassimp-dev liboctomap-dev

Build
-----

From the repository root:

.. code-block:: bash

	cmake -S . -B build
	cmake --build build -j$(nproc)
	cd build && ctest

The conic solver backends are optional and off by default:

.. code-block:: bash

	cmake -S . -B build -DCARDILLO_WITH_QOCO=ON -DCARDILLO_WITH_CLARABEL=ON -DCARDILLO_WITH_CONICXX=ON -DCARDILLO_WITH_SCS=ON

.. list-table:: CMake options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Default
     - Meaning
   * - ``CARDILLO_WITH_CLARABEL``, ``_QOCO``, ``_CONICXX``, ``_SCS``
     - OFF
     - Optional conic solver backends (fetched and built)
   * - ``CARDILLO_WITH_MOSEK``
     - OFF
     - MOSEK backend; needs a MOSEK installation and license (below)
   * - ``CARDILLO_BUILD_EXAMPLES``, ``CARDILLO_BUILD_TESTS``
     - ON when top level
     - Scene executables, tests
   * - ``CARDILLO_INSTALL``
     - ON when top level
     - Install rules and the CMake package
   * - ``CARDILLO_ENABLE_NATIVE_ARCH``
     - ON when top level
     - ``-march=native``; propagated to code using cardillo
   * - ``CARDILLO_USE_SYSTEM_COAL``
     - OFF
     - Use an installed Coal built from the fork instead of fetching it
   * - ``CARDILLO_ENABLE_TSAN``, ``CARDILLO_ENABLE_ASAN``
     - OFF
     - ThreadSanitizer, AddressSanitizer + UndefinedBehaviorSanitizer

MOSEK is commercial and not fetched. Install it and its license, then configure with
``-DCARDILLO_WITH_MOSEK=ON``; ``~/mosek/<version>`` is found automatically, other locations with
``-DMOSEK_ROOT=<dir>`` (or the environment variable ``MOSEK_HOME``):

.. code-block:: bash

	curl -LO https://download.mosek.com/stable/11.2.5/mosektoolslinux64x86.tar.bz2
	tar -xjf mosektoolslinux64x86.tar.bz2 -C ~
	cp /path/to/mosek.lic ~/mosek/mosek.lic        # or set MOSEKLM_LICENSE_FILE

Run an example
--------------

Each example scene is driven by a config file that sets the output folder, solver and time step.
Run the scenes from the repository root:

.. code-block:: bash

	build/examples/wilberforce examples/scenes/wilberforce/scene.config

If you want a compact code walkthrough after building the project, see the
worked pendulum example in :doc:`using_the_engine/example`.

Use CardilloCxx in another project
----------------------------------

Installed:

.. code-block:: bash

	cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/opt/cardillo
	cmake --build build -j$(nproc) && cmake --install build

.. code-block:: cmake

	find_package(cardillo REQUIRED)  # with CMAKE_PREFIX_PATH=$HOME/opt/cardillo
	target_link_libraries(app PRIVATE cardillo::cardillo)

The fetched dependencies (Eigen, EnTT, Coal, enabled conic backends) are installed into the same
prefix. As a subproject, examples, tests and install rules are off by default:

.. code-block:: cmake

	include(FetchContent)
	FetchContent_Declare(cardillo GIT_REPOSITORY https://github.com/CardilloCxx/CardilloCxx.git GIT_TAG main)
	FetchContent_MakeAvailable(cardillo)
	target_link_libraries(app PRIVATE cardillo::cardillo)

With ``CARDILLO_ENABLE_NATIVE_ARCH`` the library is compiled with ``-march=native``, and code using
it gets the same flag: the alignment of fixed-size Eigen members differs between instruction sets and
is part of the ABI of the headers.

Docs build
----------

The documentation is generated with Sphinx, Breathe, and Exhale. Rebuild it from ``docs/`` after changing the public headers or chapter sources:

.. code-block:: bash

	uv run make html

For details on all types (vectors, matrices, shapes, components), see :doc:`types`.
