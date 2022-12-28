#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_OPS_H
#define TINYMESH_OPS_H

#include "core/api.h"
#include "core/mesh.h"
#include "core/utils.h"

#define EIGEN_ENABLE_SPARSE
#include "core/eigen.h"

namespace tinymesh {

/**
 * Hausdorff distance of two meshes 
 */
TINYMESH_API double getHausdorffDistance(const Mesh &m0, const Mesh &m1);

/**
 * Compute per-vertex shortest distances 
 */
TINYMESH_API EigenVector getPerVertexShortestDistances(const Mesh &src, const Mesh &dst);

//! Enumerator for mesh Laplacian types
enum MeshLaplace : int {
    Adjacent,
    Cotangent,
    Belkin08,
};

/**
 * 
 */
TINYMESH_API EigenSparseMatrix getMeshLaplacian(const Mesh &mesh, MeshLaplace type);

/**
 * 
 */
TINYMESH_API EigenMatrix getHeatKernelSignatures(const EigenSparseMatrix &L, int K = 300, int nTimes = 100);

/**
 * Compute per-vertex curvature tensor (i.e., the 2nd fundamental form)
 */
TINYMESH_API void getCurvatureTensors(const Mesh &mesh, std::vector<EigenMatrix2> &tensors,
                                      const std::vector<LocalFrame> &frames);

/**
 * Compute principal curvatures and their directions [Rusinkiewicz 2004]
 */
TINYMESH_API std::tuple<EigenVector, EigenVector, EigenMatrix, EigenMatrix> getPrincipalCurvatures(
    const Mesh &mesh, bool smoothTensors = false);

/**
 * Compute principal curvatures with derivatives [Rusinkiewicz 2004] 
 */
TINYMESH_API std::tuple<EigenVector, EigenVector, EigenVector, EigenVector, EigenMatrix, EigenMatrix>
getPrincipalCurvaturesWithDerivatives(const Mesh &mesh, bool smoothTensors = false);

/**
 * Compute feature line field [Leifman and Tal 2012]
 */
TINYMESH_API EigenMatrix getFeatureLineField(const Mesh &mesh, bool smoothTensors = false);

}  // namespace tinymesh

#endif  // TINYMESH_OPS_H
