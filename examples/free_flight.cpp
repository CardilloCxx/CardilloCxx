#include "cardillo.hpp"
#include <mpi.h>
#include <petsc.h>

using namespace cardillo;

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

    real_t mass = 3.14;
    PointMass point_mass = PointMass(mass);

    

    // TODO:
    // - how to create a system of different particles
    // - distribut indices for each particle with start and offset (see misc/dofs.hpp)
    // - gather and scatter function on system level or global vectors for q and u and bodies have views insides these vectors?

    PetscFinalize();
    MPI_Finalize();
    return 0;
}
