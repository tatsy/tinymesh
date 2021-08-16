#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_REMESH_H
#define TINYMESH_REMESH_H

#include "core/api.h"

namespace tinymesh {

/**
 * Triangular remeshing [Hoppe 1996]
 */
TINYMESH_API void remeshTriangular(Mesh &mesh, double shortLength = 0.8, double longLength = 1.333, double angleThresh = 0.2, int maxiter = 5);

/**
 * QEM-based simplification [Garland and Heckbert 1997]
 */
TINYMESH_API void simplifyQEM(Mesh &mesh, int numTarget);

/**
 * Minimum direhdral hole filling [Leipa 2003]
 */
TINYMESH_API void holeFill(Mesh &mesh, double dihedralBound = Pi);

}  // namespace tinymesh

#endif  // TINYMESH_REMESH_H
