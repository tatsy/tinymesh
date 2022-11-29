#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_RESTORE_H
#define TINYMESH_RESTORE_H

#include "core/api.h"

namespace tinymesh {

/**
 * Minimum direhdral hole filling [Leipa 2003]
 * @details The triangulation algorithm is based on the following papers.
 * Barequet and Sharir, "Filling Gaps in the Boundary of a Polyhedron", 1995.
 * Liepa, "Filling Holes in Meshes", 2003. If "dihedralBound" is specified as "Pi",
 * then this method works as Barequet's method; otherwise, it works as Liepa's method.
 */
TINYMESH_API void holeFillMinDihedral(Mesh &mesh, Face *f, double dihedralBound = Pi);

/**
 * Advancing front hole filling [Zhao et al. 2007]
 */
TINYMESH_API void holeFillAdvancingFront(Mesh &mesh, Face *f);

/**
 * Context-based Coherent Surface Completion [Harary et al. 2016]
 */
TINYMESH_API void holeFillingContextCoherent(Mesh &mesh, Face *f);

}  // namespace tinymesh

#endif  // TINYMESH_RESTORE_H
