#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_HOLE_FILL_H
#define TINYMESH_HOLE_FILL_H

#include "core/common.h"

namespace tinymesh {

TINYMESH_API void hole_fill(Mesh &mesh, double dihedralBound = Pi);

}  // namespace tinymesh

#endif  // TINYMESH_HOLE_FILL_H
