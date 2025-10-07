#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <petsc.h>
// #include <Eigen/SparseCore>
// #include <cardillo/miscellaneous/coo.hpp>
// #include <unsupported/Eigen/CXX11/Tensor>

namespace cardillo {

// TODO:
// - we have to set these values by CMake
// - these values should coincide with PetscInt and PetscReal, so we can 
//   possibly use these.
// #define REAL_TYPE DOUBLE
// #define INDEX_SIZE 32
#define REAL_TYPE PetscReal
#define INDEX_SIZE PetscInt

#if REAL_TYPE == DOUBLE
using real_t = double;
#elif REAL_TYPE == LONG_DOUBLE
using real_t = long double;
#else
using real_t = float;
#endif

#if INDEX_SIZE == 8
using index_t = int8_t;
#elif INDEX_SIZE == 16
using index_t = int16_t;
#elif INDEX_SIZE == 32
using index_t = int32_t;
#else
using index_t = int64_t;
#endif

// fixed size objects
using Vector0r = Eigen::Vector<real_t, 0>;
using Matrix30r = Eigen::Matrix<real_t, 3, 0>;
// using Tensor0r = Eigen::TensorFixedSize<real_t, Eigen::Sizes<0>>;
// using Tensor300r = Eigen::TensorFixedSize<real_t, Eigen::Sizes<3, 0, 0>>;
// using Tensor330r = Eigen::TensorFixedSize<real_t, Eigen::Sizes<3, 3, 0>>;

using Vector1r = Eigen::Matrix<real_t, 1, 1>;
using Matrix11r = Eigen::Matrix<real_t, 1, 1>;

using Array2r = Eigen::Array<real_t, 2, 1>;
using Array2i = Eigen::Array<index_t, 2, 1>;
using Vector2r = Eigen::Matrix<real_t, 2, 1>;
using Vector2i = Eigen::Matrix<index_t, 2, 1>;

using Vector3r = Eigen::Matrix<real_t, 3, 1>;
using Array3i = Eigen::Array<index_t, 3, 1>;
using Matrix33r = Eigen::Matrix<real_t, 3, 3>;
using SkewSymmetricMatrix3r = Eigen::SkewSymmetricMatrix3<real_t>;
// using Tensor333r = Eigen::TensorFixedSize<real_t, Eigen::Sizes<3, 3, 3>>;

using AngleAxis3r = Eigen::AngleAxis<real_t>;
using Quaternion4r = Eigen::Quaternion<real_t>;

template <index_t RowsAtCompileTime>
using Vectorr = Eigen::Matrix<real_t, RowsAtCompileTime, 1>;

template <index_t RowsAtCompileTime>
using Arrayri = Eigen::Array<index_t, RowsAtCompileTime, 1>;
template <index_t RowsAtCompileTime>
using Vectori = Eigen::Matrix<index_t, RowsAtCompileTime, 1>;
template <index_t RowsAtCompileTime, index_t ColsAtCompileTime>
using Matrixr = Eigen::Matrix<real_t, RowsAtCompileTime, ColsAtCompileTime>;

template <index_t RowsAtCompileTime, index_t ColsAtCompileTime>
using Matrixi = Eigen::Matrix<index_t, RowsAtCompileTime, ColsAtCompileTime>;

// dynamic size objects
using ArrayXi = Eigen::Array<index_t, Eigen::Dynamic, 1>;
using ArrayXr = Eigen::Array<real_t, Eigen::Dynamic, 1>;
using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using VectorXi = Eigen::Matrix<index_t, Eigen::Dynamic, 1>;
using VectorXr = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using ArrayXXr = Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixXXi = Eigen::Matrix<index_t, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixXXr = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>;
// using TensorXr = Eigen::Tensor<real_t, 1>;
// using TensorXXr = Eigen::Tensor<real_t, 2>;
// using TensorXXXr = Eigen::Tensor<real_t, 3>;
// using TensorXXXXr = Eigen::Tensor<real_t, 4>;
// using TensorXXXXXr = Eigen::Tensor<real_t, 5>;

// reference type, read only
// TODO: Make this
// using ConstRefVectorXr = const Eigen::Ref<const VectorXr>;
// and allow for
// using RefVectorXr = Eigen::Ref<VectorXr>;
// see
// https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class
// and https://libeigen.gitlab.io/docs/classEigen_1_1Ref.html for correct usage.
using RefVectorXr = Eigen::Ref<const VectorXr>;
using RefVectorXi = Eigen::Ref<const VectorXi>;
using RefVectorXb = Eigen::Ref<const VectorXb>;
using RefMatrixXXr = Eigen::Ref<const MatrixXXr>;
using RefMatrixXXi = Eigen::Ref<const MatrixXXi>;
// using RefTensorXr = Eigen::TensorRef<const TensorXr>;
// using RefTensorXXr = Eigen::TensorRef<const TensorXXr>;
// using RefTensorXXXr = Eigen::TensorRef<const TensorXXXr>;

// maps
using MapArrayXi = Eigen::Map<ArrayXi>;
using MapVectorXi = Eigen::Map<VectorXi>;
using MapVectorXr = Eigen::Map<VectorXr>;
using MapMatrixXXr = Eigen::Map<MatrixXXr>;
// using MapTensorXr = Eigen::TensorMap<TensorXr>;
// using MapTensorXXr = Eigen::TensorMap<TensorXXr>;
// using MapQuaternion4r = Eigen::Map<const Quaternion4r>;

// stl vectors
using StlVectorXr = std::vector<real_t>;
using StlVectorXi = std::vector<index_t>;
using StlVectorStlVectorXi = std::vector<std::vector<index_t>>;
using StlVectorVectorXr = std::vector<VectorXr>;

// sparse matrices (see https://petsc.org/release/developers/matrices/)
// TODO: Use the best suited PETSc sparse matrix formats here
// // TODO: do we want CooMatrixXXr, Tripletr, SparseMatrixXXr, CscMatrixXXr,
// // CsrMatrixXXr?
// using CooMatrix = miscellaneous::CooMatrix<real_t>;
// using Triplet = Eigen::Triplet<real_t>;
// template <int StorageLayout>
// using SparseMatrix = Eigen::SparseMatrix<real_t, StorageLayout>;
// using CscMatrix = SparseMatrix<Eigen::ColMajor>;
// using CsrMatrix = SparseMatrix<Eigen::RowMajor>;

}  // namespace cardillo