# Matrix block-ops benchmark

This micro-benchmark compares a few ways to represent and operate on a sparse block matrix like the contact Jacobian used in the physics engine. Each contact contributes two 3x3 blocks coupling two bodies.

Implementations:
- EigenDense: single dense matrix (3*nc) x (3*nb) using Eigen::MatrixX<double>
- EigenSparse: sparse matrix using Eigen::SparseMatrix<double> with triplets
- BlockArray: custom structure storing two 3x3 blocks per contact row (Bi, Bj) and the indices of the two bodies; fast, cache-friendly kernels

Operations measured (looped for N iters):
- A*x: contact-space y = A * x, where x is per-body stacked 3-vectors
- rowDot: accumulate per-row dot products with x (proxy for constraint evaluation)
- At*w (acc): y += A^T * w accumulation into body-space
- At*A*x: normal-equation multiply used in CG-like solvers

Build:
This target is integrated into the top-level CMake.

Run examples:
```
# Block-array
./build/bin/matrix_bench --impl block --nb 2000 --nc 4000 --iters 10

# Eigen sparse
./build/bin/matrix_bench --impl sparse --nb 2000 --nc 4000 --iters 10

# Eigen dense (slow for large sizes)
./build/bin/matrix_bench --impl dense --nb 2000 --nc 4000 --iters 5
```

Flags:
- --impl [block|sparse|dense]
- --nb <num bodies>
- --nc <num contacts>
- --iters <iterations of each op>
- --seed <rng seed>

Notes:
- Block size is fixed at 3 for this benchmark to match linear velocity components; easy to generalize.
- The assembled blocks are diagonal 3x3 with simple deterministic scaling; swap in your real per-contact 3x3 blocks to mirror production math more closely.
- All calculations use double to keep this benchmark independent from PETSc/cardillo typedefs.