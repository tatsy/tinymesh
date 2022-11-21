#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_EIGEN_H
#define TINYMESH_EIGEN_H

using FloatType = double;
using IndexType = int64_t;

#include <Eigen/Core>
using EigenVector = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;
using EigenMatrix = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic>;

#ifdef EIGEN_ENABLE_SPARSE
#include <Eigen/SparseCore>
using EigenTriplet = Eigen::Triplet<FloatType, IndexType>;
using EigenSparseVector = Eigen::SparseVector<FloatType>;
using EigenSparseMatrix = Eigen::SparseMatrix<FloatType>;
#endif  // EIGEN_ENABLE_SPARSE

#endif  // TINYMESH_EIGEN_H
