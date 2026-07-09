#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <vector>

#if defined _WIN32 || defined __CYGWIN__
    #define CARDILLO_API __declspec(dllexport)
#else
  #define CARDILLO_API __attribute__ ((visibility ("default")))
#endif

namespace cardillo {

#define real_t double
#define index_t int

// fixed size objects
using Vector0r = Eigen::Vector<real_t, 0>;
using Matrix30r = Eigen::Matrix<real_t, 3, 0>;

using Vector1r = Eigen::Matrix<real_t, 1, 1>;
using Matrix11r = Eigen::Matrix<real_t, 1, 1>;

using Array2r = Eigen::Array<real_t, 2, 1>;
using Array2i = Eigen::Array<index_t, 2, 1>;
using Vector2r = Eigen::Matrix<real_t, 2, 1>;
using Vector2i = Eigen::Matrix<index_t, 2, 1>;

using Vector3r = Eigen::Matrix<real_t, 3, 1>;
using Vector3b = Eigen::Matrix<bool, 3, 1>;
using Array3i = Eigen::Array<index_t, 3, 1>;
using Array3b = Eigen::Array<bool, 3, 1>;
using Matrix33r = Eigen::Matrix<real_t, 3, 3>;
using SkewSymmetricMatrix3r = Eigen::SkewSymmetricMatrix3<real_t>;

using Vector4r = Eigen::Matrix<real_t, 4, 1>;
using Matrix44r = Eigen::Matrix<real_t, 4, 4>;

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

// reference type, read only
// TODO: Make this
// using ConstRefVectorXr = const Eigen::Ref<const VectorXr>;
// and allow for using RefVectorXr = Eigen::Ref<VectorXr>;
// see https://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class
// and https://libeigen.gitlab.io/docs/classEigen_1_1Ref.html for correct usage.
using RefVectorXr = Eigen::Ref<const VectorXr>;
using RefVectorXi = Eigen::Ref<const VectorXi>;
using RefVectorXb = Eigen::Ref<const VectorXb>;
using RefMatrixXXr = Eigen::Ref<const MatrixXXr>;
using RefMatrixXXi = Eigen::Ref<const MatrixXXi>;

// maps
using MapArrayXi = Eigen::Map<ArrayXi>;
using MapVectorXi = Eigen::Map<VectorXi>;
using MapVectorXr = Eigen::Map<VectorXr>;
using MapMatrixXXr = Eigen::Map<MatrixXXr>;
// using MapQuaternion4r = Eigen::Map<const Quaternion4r>;

// sparse matrices
// // TODO: do we want CooMatrixXXr, Tripletr, SparseMatrixXXr, CscMatrixXXr,
// // CsrMatrixXXr?
// using CooMatrix = miscellaneous::CooMatrix<real_t>;
// using Triplet = Eigen::Triplet<real_t>;
template <int StorageLayout>
using SparseMatrix = Eigen::SparseMatrix<real_t, StorageLayout>;
using CscMatrix = SparseMatrix<Eigen::ColMajor>;
using CsrMatrix = SparseMatrix<Eigen::RowMajor>;

}  // namespace cardillo