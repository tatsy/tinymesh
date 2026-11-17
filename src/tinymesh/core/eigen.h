#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_EIGEN_H
#define TINYMESH_EIGEN_H

#include <cstdint>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "api.h"

namespace tinymesh {

using FloatType = double;
using IndexType = std::int64_t;

using EigenVector2 = Eigen::Matrix<FloatType, 2, 1>;
using EigenVector3 = Eigen::Matrix<FloatType, 3, 1>;
using EigenVector4 = Eigen::Matrix<FloatType, 4, 1>;
using EigenVector = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
using EigenMatrix2 = Eigen::Matrix<FloatType, 2, 2>;
using EigenMatrix3 = Eigen::Matrix<FloatType, 3, 3>;
using EigenMatrix4 = Eigen::Matrix<FloatType, 4, 4>;
using EigenMatrix = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;

using EigenTriplet = Eigen::Triplet<FloatType, IndexType>;
using EigenSparseVector = Eigen::SparseVector<FloatType>;
using EigenSparseMatrix = Eigen::SparseMatrix<FloatType>;

//! SVD using Eigen
TINYMESH_API void eigenSVD(const EigenMatrix &A, EigenMatrix &U, EigenVector &Sigma, EigenMatrix &V);

}  // namespace tinymesh

#endif  // TINYMESH_EIGEN_H
