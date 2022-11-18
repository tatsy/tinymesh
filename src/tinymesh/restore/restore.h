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
TINYMESH_API void holeFill(Mesh &mesh, double dihedralBound = Pi);

}  // namespace tinymesh

#endif  // TINYMESH_RESTORE_H
