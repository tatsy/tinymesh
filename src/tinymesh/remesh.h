#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_REMESH_H
#define TINYMESH_REMESH_H

#include "common.h"

namespace tinymesh {

TINYMESH_EXPORTS void simplify(Mesh &mesh, int maxiter = 5);

}  // namespace tinymesh

#endif  // TINYMESH_REMESH_H
