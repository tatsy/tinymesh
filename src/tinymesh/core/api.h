#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_API_H
#define TINYMESH_API_H

// -----------------------------------------------------------------------------
// Forward declaration
// -----------------------------------------------------------------------------
namespace tinymesh {

class Vertex;
class Halfedge;
class Face;
class Mesh;

}  // namespace tinymesh

// -----------------------------------------------------------------------------
// API export macro
// -----------------------------------------------------------------------------

#if !defined(TINYMESH_SHARED)
#define TINYMESH_API

#elif defined(_WIN32) || defined(__CYGWIN__)
#if defined(TINYMESH_API_EXPORT)
#define TINYMESH_API __declspec(dllexport)
#else
#define TINYMESH_API __declspec(dllimport)
#endif

#elif defined(__GNUC__) || defined(__clang__)
#define TINYMESH_API __attribute__((visibility("default")))
#else
#define TINYMESH_API
#endif

#endif
