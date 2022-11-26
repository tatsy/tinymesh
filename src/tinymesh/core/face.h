#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_FACE_H
#define TINYMESH_FACE_H

#include <memory>

#include "core/api.h"
#include "core/vec.h"

namespace tinymesh {

class TINYMESH_API Face {
public:
    // Forward declaration
    class VertexIterator;
    class FaceIterator;

public:
    Face();
    Face(const Face &face) = default;
    Face(Face &&face) noexcept = default;
    virtual ~Face() = default;

    Face &operator=(const Face &face) = default;
    Face &operator=(Face &&face) noexcept = default;

    bool operator==(const Face &other) const;

    VertexIterator v_begin();
    VertexIterator v_end();
    FaceIterator f_begin();
    FaceIterator f_end();

    int index() const {
        return index_;
    }

    Vec3 normal();
    double area();

    bool isBoundary();
    bool isStatic();

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

}  // namespace tinymesh

#endif  // TINYMESH_FACE_H
