#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_EIGEN_H
#define TINYMESH_EIGEN_H

#include <cstdint>

using FloatType = double;
using IndexType = int64_t;

#include <Eigen/Core>
using EigenVector2 = Eigen::Matrix<FloatType, 2, 1>;
using EigenVector3 = Eigen::Matrix<FloatType, 3, 1>;
using EigenVector4 = Eigen::Matrix<FloatType, 4, 1>;
using EigenVector = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
using EigenMatrix2 = Eigen::Matrix<FloatType, 2, 2>;
using EigenMatrix3 = Eigen::Matrix<FloatType, 3, 3>;
using EigenMatrix4 = Eigen::Matrix<FloatType, 4, 4>;
using EigenMatrix = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;

//! SVD using Eigen
void eigenSVD(const EigenMatrix &A, EigenMatrix &U, EigenVector &Sigma, EigenMatrix &V);

// #ifdef EIGEN_ENABLE_SPARSE
#include <Eigen/SparseCore>
using EigenTriplet = Eigen::Triplet<FloatType, IndexType>;
using EigenSparseVector = Eigen::SparseVector<FloatType>;
using EigenSparseMatrix = Eigen::SparseMatrix<FloatType>;
// #endif  // EIGEN_ENABLE_SPARSE

#endif  // TINYMESH_EIGEN_H
