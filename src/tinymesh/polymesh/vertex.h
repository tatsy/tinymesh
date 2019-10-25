#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_VERTEX_H
#define TINYMESH_VERTEX_H

#include <memory>

#include "core/common.h"
#include "core/vec.h"

namespace tinymesh {

class TINYMESH_API Vertex {
public:
    // Forward declaration
    class VertexIterator;
    class InHalfedgeIterator;
    class OutHalfedgeIterator;
    class FaceIterator;

public:
    Vertex();
    Vertex(const Vertex &v) = default;
    Vertex(Vertex &&p) noexcept = default;
    explicit Vertex(const Vec3 &v);
    virtual ~Vertex() = default;

    Vertex &operator=(const Vertex &p) = default;
    Vertex &operator=(Vertex &&p) noexcept = default;

    bool operator==(const Vertex &other) const;

    int degree();
    Vec3 normal();

    VertexIterator v_begin();
    VertexIterator v_end();
    InHalfedgeIterator ihe_begin();
    InHalfedgeIterator ihe_end();
    OutHalfedgeIterator ohe_begin();
    OutHalfedgeIterator ohe_end();
    FaceIterator f_begin();
    FaceIterator f_end();

    Vec3 pos() const { return pos_; }
    void setPos(const Vec3 &pos) { pos_ = pos; }

    int index() const { return index_; }

    bool isBoundary();

private:
    Vec3 pos_;
    Halfedge *halfedge_ = nullptr;
    int index_ = -1;

    friend class Mesh;
};

// ----------
// VertexIterator
// ----------

class Vertex::VertexIterator {
public:
    explicit VertexIterator(Halfedge *he);
    bool operator!=(const VertexIterator &it) const;
    Vertex &operator*();
    Vertex *ptr() const;
    Vertex *operator->() const;
    VertexIterator &operator++();
    VertexIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

// ----------
// InHalfedgeIterator
// ----------

class Vertex::InHalfedgeIterator {
public:
    explicit InHalfedgeIterator(Halfedge *he);
    bool operator!=(const InHalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *ptr() const;
    Halfedge *operator->() const;
    InHalfedgeIterator &operator++();
    InHalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

// ----------
// OutHalfedgeIterator
// ----------

class Vertex::OutHalfedgeIterator {
public:
    OutHalfedgeIterator(Halfedge *he);
    bool operator!=(const OutHalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *ptr() const;
    Halfedge *operator->() const;
    OutHalfedgeIterator &operator++();
    OutHalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

// ----------
// FaceIterator
// ----------

class Vertex::FaceIterator {
public:
    FaceIterator(Halfedge *he);
    bool operator!=(const FaceIterator &it) const;
    Face &operator*();
    Face *ptr() const;
    Face *operator->() const;
    FaceIterator &operator++();
    FaceIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

}  // namespace tinymesh

#endif  // TINYMESH_VERTEX_H
