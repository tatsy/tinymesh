#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_COMMON_H
#define TINYMESH_COMMON_H

#include <cmath>
#include <cstdlib>
#include <iostream>

static const double Pi = 4.0 * std::atan(1.0);

// -----------------------------------------------------------------------------
// Forward declaration
// -----------------------------------------------------------------------------
namespace tinymesh {

class Vertex;
class Edge;
class Halfedge;
class Face;
class Mesh;

}  // namespace tinymesh

// -----------------------------------------------------------------------------
// API export macro
// -----------------------------------------------------------------------------

#if (defined(WIN32) || defined(_WIN32) || defined(WINCE) || defined(__CYGWIN__))
#   if defined(TINYMESH_API_EXPORT)
#       define TINYMESH_API __declspec(dllexport)
#       define TINYMESH_IMPORTS
#   else
#       define TINYMESH_API
#       define TINYMESH_IMPORTS __declspec(dllimport)
#   endif
#elif defined(__GNUC__) && __GNUC__ >= 4
#   define TINYMESH_API __attribute__((visibility ("default")))
#   define TINYMESH_IMPORTS
#else
#   define TINYMESH_API
#   define TINYMESH_IMPORTS
#endif

// -----------------------------------------------------------------------------
// Assertion with message
// -----------------------------------------------------------------------------

#ifndef __FUNCTION_NAME__
#if defined(_WIN32) || defined(__WIN32__)
#define __FUNCTION_NAME__ __FUNCTION__
#else
#define __FUNCTION_NAME__ __func__
#endif
#endif

#ifndef DISABLE_ASSERT
#define Assertion(PREDICATE, ...) \
do { \
    if (!(PREDICATE)) { \
        std::cerr << "Asssertion \"" \
        << #PREDICATE << "\" failed in " << __FILE__ \
        << " line " << __LINE__ \
        << " in function \"" << (__FUNCTION_NAME__) << "\"" \
        << " : "; \
        fprintf(stderr, __VA_ARGS__); \
        std::cerr << std::endl; \
        std::abort(); \
    } \
} while (false)
#else  // NDEBUG
#define Assertion(PREDICATE, ...) do {} while (false)
#endif // NDEBUG



#endif // TINYMESH_COMMON_H
