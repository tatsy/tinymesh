#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_VERTEX_H
#define TINYMESH_VERTEX_H

#include <memory>

#include "common.h"
#include "vector.h"

namespace tinymesh {

class TINYMESH_EXPORTS Vertex {
public:
    // Forward declaration
    class VertexIterator;
    class InHalfedgeIterator;
    class OutHalfedgeIterator;

public:
    Vertex();
    Vertex(const Vertex &p) = default;
    Vertex(Vertex &&p) noexcept = default;
    explicit Vertex(const Vector &v);
    virtual ~Vertex() = default;

    Vertex &operator=(const Vertex &p) = default;
    Vertex &operator=(Vertex &&p) noexcept = default;

    VertexIterator v_begin();
    VertexIterator v_end();
    InHalfedgeIterator ihe_begin();
    InHalfedgeIterator ihe_end();
    OutHalfedgeIterator ohe_begin();
    OutHalfedgeIterator ohe_end();

    Vector pt() const { return m_pt; }
    void setPt(const Vector &pt) { m_pt = pt; }

private:
    Vector m_pt;
    Halfedge *m_he;
    Face *m_face;

    friend class Mesh;
};

// ----------
// VertexIterator
// ----------

class Vertex::VertexIterator {
public:
    VertexIterator(Halfedge *he);
    bool operator!=(const VertexIterator &it) const;
    Vertex &operator*();
    Vertex *operator->() const;
    VertexIterator &operator++();
    VertexIterator operator++(int);

private:
    Halfedge *m_he, *m_init;
};

// ----------
// InHalfedgeIterator
// ----------

class Vertex::InHalfedgeIterator {
public:
    InHalfedgeIterator(Halfedge *he);
    bool operator!=(const InHalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *operator->() const;
    InHalfedgeIterator &operator++();
    InHalfedgeIterator operator++(int);

private:
    Halfedge *m_he, *m_init;
};

// ----------
// OutHalfedgeIterator
// ----------

class Vertex::OutHalfedgeIterator {
public:
    OutHalfedgeIterator(Halfedge *he);
    bool operator!=(const OutHalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *operator->() const;
    OutHalfedgeIterator &operator++();
    OutHalfedgeIterator operator++(int);

private:
    Halfedge *m_he, *m_init;
};

}  // namespace tinymesh

#endif  // TINYMESH_VERTEX_H
