#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_SMOOTH_H
#define TINYMESH_SMOOTH_H

#include "polymesh/mesh.h"

namespace tinymesh {

TINYMESH_API void smooth(Mesh &mesh, double strength = 1.0);

}  // namespace tinymesh

#endif  // TINYMESH_SMOOTH_H
