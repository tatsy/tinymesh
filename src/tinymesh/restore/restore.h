#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_RESTORE_H
#define TINYMESH_RESTORE_H

#include "core/api.h"
#include "core/mesh.h"
#include "core/face.h"
#include "core/constants.h"

namespace tinymesh {

/**
 * Minimum direhdral hole filling [Liepa 2003]
 * @details The triangulation algorithm is based on the following papers:
 * Barequet and Sharir, "Filling Gaps in the Boundary of a Polyhedron", 1995.
 * Liepa, "Filling Holes in Meshes", 2003.
 * 
 * This method works as Barequet's method if "dihedralBound == kPi",
 * and as Liepa's method if dihedralBound == 0.0.
 * Otherwise, the method uses "dihedralBound" as possible dihedral angle increase.
 */
TINYMESH_API void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound = kPi);

/**
 * Advancing front hole filling [Zhao et al. 2007]
 */
TINYMESH_API void holeFillAdvancingFront(Mesh &mesh, Face *face);

/**
 * Context-based Coherent Surface Completion [Harary et al. 2016]
 */
TINYMESH_API void holeFillContextCoherent(Mesh &mesh, int maxiters = 200);

}  // namespace tinymesh

#endif  // TINYMESH_RESTORE_H
