#pragma once

#include <petscsys.h>
#include <Eigen/Dense>
#include "types.hpp"

// A simple wrapper function
namespace mylib {
    inline void say_hello(const Eigen::VectorXd& v) {
        PetscPrintf(PETSC_COMM_WORLD, "Hello from PETSc + Eigen! Vector size = %d\n", (int)v.size());
    }
}
