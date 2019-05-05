#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_SIMPLIFY_H
#define TINYMESH_SIMPLIFY_H

#include "core/common.h"

namespace tinymesh {

TINYMESH_API void simplifyIncremental(Mesh &mesh, int numTarget);

}  // namespace tinymesh

#endif  // TINYMESH_SIMPLIFY_H
