#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_REMESH_H
#define TINYMESH_REMESH_H

#include "core/common.h"

namespace tinymesh {

TINYMESH_API void remeshIncremental(Mesh &mesh, double shortLength = 0.8, double longLength = 1.333, double angleThresh = 0.2, int maxiter = 5);

}  // namespace tinymesh

#endif  // TINYMESH_REMESH_H
