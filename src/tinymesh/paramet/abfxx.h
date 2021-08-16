#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_ABFXX_H
#define TINYMESH_ABFXX_H

#include "core/common.h"

namespace tinymesh {

TINYMESH_API void abfxx(Mesh &mesh, int maxiter = 10, double epsilon = 1.0e-6);

}  // namespace tinymesh

#endif  // TINYMESH_ABFXX_H
