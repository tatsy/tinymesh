#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_UTILS_H
#define TINYMESH_UTILS_H

#include "core/api.h"
#include "core/mesh.h"

#define EIGEN_ENABLE_SPARSE
#include "core/eigen.h"

namespace tinymesh {

enum MeshLaplace : int {
    Adjacent,
    Cotangent,
    Belkin08,
};

TINYMESH_API EigenSparseMatrix getMeshLaplacian(Mesh &mesh, MeshLaplace type);

TINYMESH_API EigenMatrix getHeatKernelSignatures(Mesh &mesh, int K = 300, int nTimes = 100);

}  // namespace tinymesh

#endif  // TINYMESH_UTILS_H
