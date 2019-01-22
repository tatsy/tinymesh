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
#       define TINYMESH_EXPORTS __declspec(dllexport)
#       define TINYMESH_IMPORTS
#   else
#       define TINYMESH_EXPORTS
#       define TINYMESH_IMPORTS __declspec(dllimport)
#   endif
#elif defined(__GNUC__) && __GNUC__ >= 4
#   define TINYMESH_EXPORTS __attribute__((visibility ("default")))
#   define TINYMESH_IMPORTS
#else
#   define TINYMESH_EXPORTS
#   define TINYMESH_IMPORTS
#endif

// -----------------------------------------------------------------------------
// Message handlers
// -----------------------------------------------------------------------------

#define Info(...) \
do { \
    std::cout << "[INFO] "; \
    fprintf(stdout, __VA_ARGS__); \
    std::cerr << std::endl; \
} while (false);

#define Warning(...) \
do { \
    std::cerr << "[WARNING] "; \
    fprintf(stdout, __VA_ARGS__); \
    std::cerr << std::endl; \
} while (false);

#define FatalError(...) \
do { \
    std::cerr << "[ERROR] "; \
    fprintf(stderr, __VA_ARGS__); \
    std::cerr << std::endl; \
    std::abort(); \
} while (false);

#ifndef NDEBUG
#define Debug(...) \
do { \
    std::cerr << "[DEBUG] "; \
    fprintf(stdout, __VA_ARGS__); \
    std::cerr << std::endl; \
} while (false);
#else
#define Debug(...)
#endif

#endif // TINYMESH_COMMON_H
