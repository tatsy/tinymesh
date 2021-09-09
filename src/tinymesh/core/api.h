#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_API_H
#define TINYMESH_API_H

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

#endif // TINYMESH_API_H
