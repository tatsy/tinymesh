#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_OPENMP_H
#define TINYMESH_OPENMP_H

#if defined(_OPENMP)

#if defined(_MSC_VER)
#define TINYMESH_OMP_PARALLEL_FOR __pragma(omp parallel for) for
#else
#define TINYMESH_OMP_PARALLEL_FOR _Pragma("omp parallel for") for
#endif

#else
#define TINYMESH_OMP_PARALLEL_FOR for
#endif

#endif  // TINYMESH_OPENMP_H
