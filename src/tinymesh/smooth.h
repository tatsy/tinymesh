#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_SMOOTH_H
#define TINYMESH_SMOOTH_H

#include "mesh.h"

namespace tinymesh {

TINYMESH_EXPORTS void smooth(Mesh &mesh, double strength = 0.5);

}  // namespace tinymesh

#endif  // TINYMESH_SMOOTH_H
