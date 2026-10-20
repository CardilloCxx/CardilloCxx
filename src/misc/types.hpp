#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>

#if defined _WIN32 || defined __CYGWIN__
    #define CARDILLO_API __declspec(dllexport)
#else
  #define CARDILLO_API __attribute__ ((visibility ("default")))
#endif

namespace cardillo {

// Scalar and index types. real_t must equal Coal's scalar type (double unless Coal is built with
// COAL_USE_FLOAT_PRECISION), which is checked in collision/coal_types.hpp; the conic solvers
// (Clarabel, QOCO, ConicXX) work in double precision as well.
using real_t = double;
using index_t = int;

// fixed size objects
using Vector2r = Eigen::Matrix<real_t, 2, 1>;
using Vector3r = Eigen::Matrix<real_t, 3, 1>;
using Vector4r = Eigen::Matrix<real_t, 4, 1>;
using Vector6r = Eigen::Matrix<real_t, 6, 1>;
using Matrix33r = Eigen::Matrix<real_t, 3, 3>;
using Matrix44r = Eigen::Matrix<real_t, 4, 4>;
using Matrix66r = Eigen::Matrix<real_t, 6, 6>;
using SkewSymmetricMatrix3r = Eigen::SkewSymmetricMatrix3<real_t>;

using AngleAxis3r = Eigen::AngleAxis<real_t>;
using Quaternion4r = Eigen::Quaternion<real_t>;

template <index_t RowsAtCompileTime, index_t ColsAtCompileTime>
using Matrixr = Eigen::Matrix<real_t, RowsAtCompileTime, ColsAtCompileTime>;

// dynamic size objects
using VectorXr = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using MatrixXXr = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;

// sparse matrices
template <int StorageLayout>
using SparseMatrix = Eigen::SparseMatrix<real_t, StorageLayout>;
using CscMatrix = SparseMatrix<Eigen::ColMajor>;

}  // namespace cardillo