#ifdef _MSC_VER
#pragma once
#endif

#ifndef COMMON_OPENMP_H
#define COMMON_OPENMP_H

#if defined(_OPENMP)
#include <omp.h>
#if defined(_MSC_VER)
#define omp_pragma __pragma(omp parallel for)
#define omp_critical __pragma(omp critical)
#else
#define omp_pragma _Pragma("omp parallel for")
#define omp_critical _Pragma("omp parallel for")
#endif
#define omp_parallel_for omp_pragma for
#define omp_lock_t omp_lock_t
#define omp_init_lock(lock) omp_init_lock(lock)
#define omp_destroy_lock(lock) omp_destroy_lock(lock)
#define omp_lock(lock) omp_set_lock(lock)
#define omp_unlock(lock) omp_unset_lock(lock)
#else
#define omp_set_num_threads(n)
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_parallel_for for
#define omp_critical
#define omp_lock_t int
#define omp_init_lock(lock)
#define omp_destroy_lock(lock)
#define omp_lock(lock)
#define omp_unlock(lock)
#endif

#endif  // COMMON_OPENMP_H