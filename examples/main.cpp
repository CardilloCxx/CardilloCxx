#include "cardillo.hpp"
#include <mpi.h>
#include <petsc.h>   // full PETSc API

// run with mpirun -np 2 ./s/example 
int main(int argc, char** argv) {
    // Initialize MPI and PETSc
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    int rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if (rank == 0) {
        PetscPrintf(PETSC_COMM_SELF, 
            "Running with %d MPI processes.\n", size);
    }

    // Test Eigen + PETSc wrapper
    Eigen::VectorXd x = Eigen::VectorXd::Ones(5);
    cardillo::say_hello(x);

    PetscFinalize();
    MPI_Finalize();
    return 0;
}
