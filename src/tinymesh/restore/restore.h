#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_RESTORE_H
#define TINYMESH_RESTORE_H

#include "core/api.h"

namespace tinymesh {

/**
 * Minimum direhdral hole filling [Leipa 2003]
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
