#ifdef _MSC_VER
#pragma once
#endif 

#ifndef TINYMESH_FACE_H
#define TINYMESH_FACE_H

#include <memory>

#include "core/common.h"

namespace tinymesh {

class TINYMESH_API Face {
public:
    // Forward declaration
    class VertexIterator;

public:
    Face();
    Face(const Face &face) = default;
    Face(Face &&face) noexcept = default;
    virtual ~Face() = default;

    Face &operator=(const Face &face) = default;
    Face &operator=(Face &&face) noexcept = default;

    VertexIterator v_begin();
    VertexIterator v_end();

private:
    Halfedge *m_he;
    int m_index;

    friend class Mesh;
};

class TINYMESH_API Face::VertexIterator {
public:
    explicit VertexIterator(Halfedge *he);
    bool operator!=(const VertexIterator &it) const;
    Vertex &operator*();
    Vertex *operator->() const;
    Vertex *ptr() const;
    VertexIterator &operator++();
    VertexIterator operator++(int);

private:
    Halfedge *m_he, *m_init;
};

}  // namespace tinymesh

#endif  // TINYMESH_FACE_H
