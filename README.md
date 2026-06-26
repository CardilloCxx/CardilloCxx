# CardilloCxx

CardilloCxx is a C++ physics engine with contact handling, multiple solver backends and many example scenes.

## What this project builds

- Core static library: `cardillo`
- Main executable: `build/bin/main`
- Scene/example and benchmark targets
- Optional interior-point solvers:
   - QOCO (CPU and optional CUDA backend)
   - Clarabel (built through Clarabel.rs, requires Rust/Cargo)

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

## Configure and build

From repository root:

```bash
cd CardilloMPI
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

## Running examples

```bash
./build/bin/main ./examples/scenes/wilberforce/scene.config
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

## Build Instructions

Run the build commands from the `docs` directory:

```bash
cd docs
uv run make html
```

After the build completes, open `docs/_build/html/index.html` in your web browser to preview the site.
