#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_REMESH_H
#define TINYMESH_REMESH_H

#include "core/common.h"

namespace tinymesh {

TINYMESH_API void remeshIncremental(Mesh &mesh, double ratioLower = 0.8, double ratioUpper = 1.333, int maxiter = 5);

}  // namespace tinymesh

#endif  // TINYMESH_REMESH_H
