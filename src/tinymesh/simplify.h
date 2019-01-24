#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_SIMPLIFY_H
#define TINYMESH_SIMPLIFY_H

#include "common.h"

namespace tinymesh {

TINYMESH_EXPORTS void simplify(Mesh &mesh, double ratio = 0.5, int num_remains = -1);

}  // namespace tinymesh

#endif  // TINYMESH_SIMPLIFY_H
