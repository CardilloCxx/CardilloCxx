#pragma once

#include <petscsys.h>
#include <Eigen/Dense>
#include "misc/types.hpp"
#include "misc/dofs.hpp"
#include "physics/physics_system.hpp"

// A simple wrapper function
namespace cardillo {
    inline void say_hello(const Eigen::VectorXd& v) {
        PetscPrintf(PETSC_COMM_WORLD, "Hello from PETSc + Eigen! Vector size = %d\n", (int)v.size());
    }
}
