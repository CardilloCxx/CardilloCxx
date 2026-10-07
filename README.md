# CardilloCxx

CardilloCxx is a C++ physics engine with contact handling, multiple solver backends and many example scenes.

<p align="center">
  <img src="https://raw.githubusercontent.com/CardilloCxx/CardilloCxx/main/docs/logo.png" alt="logo" width="400"/>
</p>

## What this project builds

- Core static library: `cardillo`
- One executable per scene under `build/bin/`
- Scene/example and benchmark targets
- Optional interior-point solvers:
   - QOCO (CPU and optional CUDA backend)
   - Clarabel (built through Clarabel.rs, requires Rust/Cargo)
   - ConicXX (plain CMake/C++17 target, no Rust/CUDA required; the only
     backend that warm-starts and reuses its KKT factorization across steps
     when the active contact set is unchanged)

## Prerequisites

Required on Linux:

- CMake 3.16+
- C++17 compiler (GCC/Clang)
- Python 3

Typical Ubuntu packages:

```bash
sudo apt update
sudo apt install -y build-essential cmake ninja-build git curl pkg-config \
   python3 python3-dev python3-pip libeigen3-dev
```

Mac packages:

```bash
brew update
xcode-select --install
brew install cmake ninja git curl pkg-config python3 eigen
```

Coal also needs Assimp and Octomap headers on many Linux setups:

```bash
sudo apt update
sudo apt install libassimp-dev liboctomap-dev
```

```Bash
brew install assimp octomap
```

## Install Rust/Cargo (needed for Clarabel)

Reason: Clarabel is built via `Clarabel.rs` in the CMake build. Without Cargo, Clarabel targets cannot be compiled.

```bash
curl https://sh.rustup.rs -sSf | sh -s -- -y
source "$HOME/.cargo/env"
cargo --version
```

Optional: persist Cargo on PATH for new shells:

```bash
echo 'source "$HOME/.cargo/env"' >> ~/.bashrc
```

## CUDA dependency for QOCO CUDA backend

QOCO CUDA backend requires NVIDIA cuDSS in addition to CUDA Toolkit.

- cuDSS: https://developer.nvidia.com/cudss

If cuDSS is not found, build will continue with CPU backend only.
This is the default case for MacOS.

## Configure and build

From repository root:

```bash
cd CardilloCxx
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

## QOCO backend options

Enable/disable CUDA variant at configure time:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DQOCO_USE_CUDA_BACKEND=ON
```

Runtime backend selection is controlled by config key `qoco.backend`:

- `auto` (prefer CUDA if available)
- `cpu`
- `cuda`

## Optional MOSEK backend

`solver.type = mosek` uses the commercial [MOSEK](https://www.mosek.com) conic optimizer. It is off by
default and MOSEK is not shipped with this repository; install it and its license yourself:

```bash
# MOSEK 11 (Optimizer API) into ~/mosek/<version>
curl -LO https://download.mosek.com/stable/11.2.5/mosektoolslinux64x86.tar.bz2
tar -xjf mosektoolslinux64x86.tar.bz2 -C ~
cp /path/to/mosek.lic ~/mosek/mosek.lic        # or set MOSEKLM_LICENSE_FILE

cmake -S . -B build -DCARDILLO_WITH_MOSEK=ON   # add -DMOSEK_ROOT=<dir> for non-default locations
```

`~/mosek/<version>` is found automatically; otherwise pass `-DMOSEK_ROOT` (or set `MOSEK_HOME`) to the
installation or its `tools/platform/<platform>` directory. Tests and the solver benchmark are built
with `-DCARDILLO_BUILD_TESTS=ON` (`ctest --test-dir build`, `./build/tests/mosek_benchmark`).

## Running examples

```bash
./build/bin/wilberforce ./examples/scenes/wilberforce/scene.config
```

## Troubleshooting

Clean rebuild

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

# Documentation

The documentation is built using Sphinx and Doxygen.

## Prerequisites

Install the required system tools:

```bash
sudo apt update && sudo apt install -y doxygen graphviz
```

```bash
brew install doxygen graphviz
```

## Build Instructions

Run the build commands from the `docs` directory:

```bash
cd docs
uv run make html
```

After the build completes, open `docs/_build/html/index.html` in your web browser to preview the site.
