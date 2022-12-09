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

enum MeshLaplace : int {
    Adjacent,
    Cotangent,
    Belkin08,
};

TINYMESH_API EigenSparseMatrix getMeshLaplacian(const Mesh &mesh, MeshLaplace type);

TINYMESH_API EigenMatrix getHeatKernelSignatures(const EigenSparseMatrix &L, int K = 300, int nTimes = 100);

//! Compute per-vertex curvature tensor (i.e., the 2nd fundamental form)
TINYMESH_API void getCurvatureTensors(const Mesh &mesh, std::vector<EigenMatrix2> &tensors,
                                      const std::vector<LocalFrame> &frames);

//! Compute principal curvatures and their directions
TINYMESH_API std::tuple<std::vector<double>, std::vector<double>, std::vector<Vec3>, std::vector<Vec3>>
getPrincipalCurvatures(const Mesh &mesh);

}  // namespace tinymesh

#endif  // TINYMESH_OPS_H
