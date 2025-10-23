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
☑ Check why we need to deduplicate contacts in naive_partitioner.cpp.

☑ Implement friction.

☑ Clean up Collision detection code, make it modular.

☑ Warm start the percussive reactions by try to match them to the previous frame impulses.

☐ Try using PETsc solver and matrices.

☐ Use a bounding volume hierachy for 
partitioning and broadphase collision detection.

☐ Implement a RigidBody from a mesh.

☐ Implement sleeping for rigid bodies?

☐ Adaptive time steps?

☐ Fix MPI desync.

☐ Add output folder to config.

☐ Setup environment for dry testing convergence of different solvers.

☐ Enable loading from checkpoints.


## References
https://www.researchgate.net/profile/Dinesh-Manocha/publication/2807460_Collision_Queries_using_Oriented_Bounding_Boxes/links/56cb392008ae5488f0daea80/Collision-Queries-using-Oriented-Bounding-Boxes.pdf