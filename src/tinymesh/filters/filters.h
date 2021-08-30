#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_FILTERS_H
#define TINYMESH_FILTERS_H

#include "core/api.h"
#include "core/mesh.h"

namespace tinymesh {

/**
 * Laplacian smoothing
 */
TINYMESH_API void smoothLaplacian(Mesh &mesh, double strength = 1.0, bool cotangent_weight = false, int iterations = 3);

/**
 * Taubin smoothing [Taubin et al. 1995]
 */
TINYMESH_API void smoothTaubin(Mesh &mesh, double shrink = 0.5, double inflate = 0.53, int iterations = 3);

/**
 * Implicit fairing [Gesbrun et al. 1999]
 */
TINYMESH_API void implicitFairing(Mesh &mesh, double lambda = 1.0, int iterations = 1);

/**
 * Denoising by normal Gaussian filter [Ohtake et al. 2001]
 */
TINYMESH_API void denoiseNormalGaussian(Mesh &mesh, double sigma = 0.2, int iterations = 5);

/**
 * Denoising by normal bilateral filter [Zhen et al. 2011]
 */
TINYMESH_API void denoiseNormalBilateral(Mesh &mesh, double sigmaCenter = 0.2, double sigmaNormal = 0.1, int iterations = 5);

/**
 * Denoising by L0 smoothing [He and Schaefer 2013]
 */
TINYMESH_API void denoiseL0Smooth(Mesh &mesh, double alpha = 0.1, double beta = 1.0e-3);

}  // namespace tinymesh

#endif  // TINYMESH_FILTERS_H
