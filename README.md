# CardilloCxx

CardilloCxx is a C++20 nonsmooth multibody engine: rigid bodies, point masses, Cosserat rods and soft
bodies with frictional contact, integrated by a Moreau-θ time-stepping scheme with several solver backends.

<p align="center">
  <img src="https://raw.githubusercontent.com/CardilloCxx/CardilloCxx/main/docs/logo.png" alt="logo" width="400"/>
</p>

The build produces the shared library `libcardillo` (headers included as `<cardillo/...>`), one
executable per example scene and the tests. The conic solver backends QOCO, Clarabel, ConicXX, SCS and
MOSEK are optional.

## Prerequisites

- CMake ≥ 3.24, a C++20 compiler (GCC or Clang), git
- VTK (a system VTK is used if found, otherwise it is built from source; on Debian/Ubuntu a
  ParaView package that ships its own VTK also works and conflicts with `libvtk9-dev`)
- Boost, Assimp and Octomap development packages (needed by Coal)
- OpenMP (optional)
- Rust/cargo, only for the Clarabel backend; a MOSEK installation and license, only for MOSEK

Eigen 5 and EnTT are used from an installation if CMake finds them, otherwise they are fetched. Coal
is fetched from the fork [41pha1/coal](https://github.com/41pha1/coal), which adds the custom shapes
of the beam hull colliders. On Ubuntu:

```bash
sudo apt install build-essential cmake git libvtk9-dev libboost-dev libassimp-dev liboctomap-dev
curl https://sh.rustup.rs -sSf | sh   # only for Clarabel
```

## Build

```bash
cmake -S . -B build                    # default build type RelWithDebInfo
cmake --build build -j$(nproc)
cd build && ctest
```

The conic backends are `OFF` by default; enable them with their CMake options, e.g.

```bash
cmake -S . -B build -DCARDILLO_WITH_CONICXX=ON -DCARDILLO_WITH_QOCO=ON -DCARDILLO_WITH_CLARABEL=ON -DCARDILLO_WITH_SCS=ON
```

Selecting a backend in a scene config that was not built throws a runtime error naming the option.

| Option | Default | Meaning |
| --- | --- | --- |
| `CARDILLO_WITH_CLARABEL` | OFF | Clarabel interior-point backend (needs Rust/cargo) |
| `CARDILLO_WITH_QOCO` | OFF | QOCO interior-point backend |
| `CARDILLO_WITH_CONICXX` | OFF | ConicXX interior-point backend |
| `CARDILLO_WITH_SCS` | OFF | SCS ADMM backend |
| `CARDILLO_WITH_MOSEK` | OFF | MOSEK interior-point backend (needs MOSEK and a license, see below) |
| `CARDILLO_SCS_OPENMP` | ON | Build SCS with OpenMP-parallel cone projections |
| `CARDILLO_BUILD_EXAMPLES` | ON if top level | Build the scene executables |
| `CARDILLO_BUILD_TESTS` | ON if top level | Build the tests |
| `CARDILLO_INSTALL` | ON if top level | Install rules and the CMake package `cardillo` |
| `CARDILLO_ENABLE_NATIVE_ARCH` | ON if top level | `-march=native`, also applied to code using cardillo |
| `CARDILLO_USE_SYSTEM_COAL` | OFF | Use an installed Coal built from the fork |
| `CARDILLO_ENABLE_TSAN` | OFF | ThreadSanitizer (use with `CMAKE_BUILD_TYPE=Debug`) |
| `CARDILLO_ENABLE_ASAN` | OFF | AddressSanitizer and UndefinedBehaviorSanitizer |

### MOSEK

MOSEK is commercial and not fetched; install it and its license yourself:

```bash
curl -LO https://download.mosek.com/stable/11.2.5/mosektoolslinux64x86.tar.bz2
tar -xjf mosektoolslinux64x86.tar.bz2 -C ~
cp /path/to/mosek.lic ~/mosek/mosek.lic        # or set MOSEKLM_LICENSE_FILE
cmake -S . -B build -DCARDILLO_WITH_MOSEK=ON   # add -DMOSEK_ROOT=<dir> for other locations
```

With Clarabel and QOCO also enabled, `mosek_tests` (ctest) and `build/tests/mosek_benchmark` compare
MOSEK on cone programs against them.

## Run a scene

Run scenes from the repository root (mesh paths are relative to it):

```bash
build/examples/wilberforce examples/scenes/wilberforce/scene.config
```

Target names are listed in `examples/CMakeLists.txt`. Output is written as VTK files and a CSV of
tracked bodies, as configured in the scene's config file.

## Use CardilloCxx in another project

Installed package:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/opt/cardillo
cmake --build build -j$(nproc) && cmake --install build
```

```cmake
find_package(cardillo REQUIRED)            # CMAKE_PREFIX_PATH=$HOME/opt/cardillo
target_link_libraries(app PRIVATE cardillo::cardillo)
```

The fetched dependencies (Eigen, EnTT, Coal, enabled conic backends) are installed into the same prefix.
As a subproject (examples, tests and install rules are then off by default):

```cmake
include(FetchContent)
FetchContent_Declare(cardillo GIT_REPOSITORY https://github.com/CardilloCxx/CardilloCxx.git GIT_TAG main)
FetchContent_MakeAvailable(cardillo)
target_link_libraries(app PRIVATE cardillo::cardillo)
```

```cpp
#include <cardillo/physics/api/physics.hpp>
```

## Documentation

The documentation is built with Sphinx and Doxygen (needs `doxygen`, `graphviz` and `uv`):

```bash
cd docs && uv run make html
```

Open `docs/_build/html/index.html` to view it.
