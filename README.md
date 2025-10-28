# CardilloMPI

Minimal example how to create a library using CMake with PETSc and Eigen3 dependencies.

## Compile

```bash
mkdir build
cd build
cmake ..
make
```

## Run

```bash
mpirun -np 4 ./examples/example
```

## Benchmark
![image](benchmarks/analysis/benchmark_results.png)
![image](benchmarks/matrix_bench/bench_plot_big.png)


## Todos

☑ Implement friction.

☑ Clean up Collision detection code, make it modular.

☑ Warm start the percussive reactions by trying to match them to the previous frame impulses.

☑ Add output folder to config.

☑ Change Physics System mass matrix from MatrixXXr to VectorXr diagonal storage.

☑ Deduplicate close contacts for meshes.

☑ Binary vtk output instead of ascii for smaller file sizes.

☑ Implement a RigidBody from a mesh.

☐ Try using PETsc solver and matrices.

☐ Use a bounding volume hierarchy for
partitioning and broadphase collision detection.

☐ Implement sleeping for rigid bodies?

☐ Adaptive time steps?

☐ Fix MPI desync.

☐ Setup environment for dry testing convergence of different solvers.

☐ Enable loading from checkpoints.

☐ numBodies() can be implemented more efficiently by caching the count of bodies and updating it on addition/removal.

☐ Fix mesh scaling issue.

☐ Fix mesh collision normal direction issue. Temporarily fixed by pointing normals away from mesh centroid.

☐ Add a scene interface that sets up different scenarios and is selected from config and loaded in main.

☐ Check if instancing multiple PhysicsSystem mesh assets works efficiently.

☐ Normalize meshes sizes on load.

## References
https://www.researchgate.net/profile/Dinesh-Manocha/publication/2807460_Collision_Queries_using_Oriented_Bounding_Boxes/links/56cb392008ae5488f0daea80/Collision-Queries-using-Oriented-Bounding-Boxes.pdf

@misc{coalweb,
   author = {Jia Pan and Sachin Chitta and Dinesh Manocha and Florent Lamiraux and Joseph Mirabel and Justin Carpentier and Louis Montaut and others},
   title = {Coal: an extension of the Flexible Collision Library},
   howpublished = {https://github.com/coal-library/coal},
   year = {2015--2024}
}

## Building COAL locally

1. Clone and build COAL:

```bash
wget https://github.com/coal-library/coal/releases/download/v3.0.2/coal-3.0.2.tar.gz
tar -xvf coal-3.0.2.tar.gz
cd coal-3.0.2
curl -fsSL https://pixi.sh/install.sh | sh
pixi run test
```

This will create `lib/` and `include/` in the COAL build directory.

2. Point CardilloMPI to the COAL build when configuring CMake:

```bash
cd CardilloMPI
mkdir -p build && cd build
cmake .. -DCMAKE_PREFIX_PATH=<....>/coal-3.0.2/build
cmake --build . -- -j$(nproc)
```

This README's COAL instructions assume the repository contains a copy of COAL under `lib/coal-3.0.2`. Adjust paths if you placed it elsewhere.

