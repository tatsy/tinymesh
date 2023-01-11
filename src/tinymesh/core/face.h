#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_FACE_H
#define TINYMESH_FACE_H

#include <memory>

#include "core/api.h"
#include "core/vec.h"

namespace tinymesh {

class Triangle;

class TINYMESH_API Face {
public:
    // Forward declaration
    class VertexIterator;
    class ConstVertexIterator;
    class HalfedgeIterator;
    class ConstHalfedgeIterator;
    class FaceIterator;
    class ConstFaceIterator;

public:
    Face();
    Face(const Face &face) = default;
    Face(Face &&face) noexcept = default;
    virtual ~Face() = default;

    Face &operator=(const Face &face) = default;
    Face &operator=(Face &&face) noexcept = default;

    bool operator==(const Face &other) const;

    Triangle toTriangle() const;

    VertexIterator v_begin();
    VertexIterator v_end();
    ConstVertexIterator v_begin() const;
    ConstVertexIterator v_end() const;

    HalfedgeIterator he_begin();
    HalfedgeIterator he_end();
    ConstHalfedgeIterator he_begin() const;
    ConstHalfedgeIterator he_end() const;

    FaceIterator f_begin();
    FaceIterator f_end();
    ConstFaceIterator f_begin() const;
    ConstFaceIterator f_end() const;

    int index() const {
        return index_;
    }

    Vec3 normal() const;
    double area() const;
    int numCorners() const;

    bool isHole() const;
    bool isBoundary() const;
    bool isLocked() const;

private:
    Halfedge *halfedge_ = nullptr;
    int index_ = -1;

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
    Halfedge *halfedge_, *iter_;
};

class TINYMESH_API Face::ConstVertexIterator {
public:
    explicit ConstVertexIterator(Halfedge *he);
    bool operator!=(const ConstVertexIterator &it) const;
    const Vertex &operator*() const;
    const Vertex *operator->() const;
    const Vertex *ptr() const;
    ConstVertexIterator &operator++();
    ConstVertexIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

class TINYMESH_API Face::HalfedgeIterator {
public:
    explicit HalfedgeIterator(Halfedge *he);
    bool operator!=(const HalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *operator->() const;
    Halfedge *ptr() const;
    HalfedgeIterator &operator++();
    HalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

class TINYMESH_API Face::ConstHalfedgeIterator {
public:
    explicit ConstHalfedgeIterator(Halfedge *he);
    bool operator!=(const ConstHalfedgeIterator &it) const;
    const Halfedge &operator*() const;
    const Halfedge *operator->() const;
    const Halfedge *ptr() const;
    ConstHalfedgeIterator &operator++();
    ConstHalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

class TINYMESH_API Face::FaceIterator {
public:
    explicit FaceIterator(Halfedge *he);
    bool operator!=(const FaceIterator &it) const;
    Face &operator*();
    Face *operator->() const;
    Face *ptr() const;
    FaceIterator &operator++();
    FaceIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

class TINYMESH_API Face::ConstFaceIterator {
public:
    explicit ConstFaceIterator(Halfedge *he);
    bool operator!=(const ConstFaceIterator &it) const;
    const Face &operator*() const;
    const Face *operator->() const;
    const Face *ptr() const;
    ConstFaceIterator &operator++();
    ConstFaceIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

}  // namespace tinymesh

#endif  // TINYMESH_FACE_H
