#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_SMOOTH_H
#define TINYMESH_SMOOTH_H

#include "polymesh/mesh.h"

namespace tinymesh {

//! Laplacian smoothing
TINYMESH_API void laplace_smooth(Mesh &mesh, double strength = 1.0, int iterations = 3);

//! Taubin smoothing [Taubin et al. 1995]
TINYMESH_API void taubin_smooth(Mesh &mesh, double shrink = 1.0, double inflate = 1.0, int iterations = 3);

//! Implicit fairing [Gesbrun et al. 1999]
TINYMESH_API void implicit_fair(Mesh &mesh, double lambda = 1.0, int iterations = 1);

}  // namespace tinymesh

#endif  // TINYMESH_SMOOTH_H
