#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_VERTEX_H
#define TINYMESH_VERTEX_H

#include <memory>

#include "core/api.h"
#include "core/vec.h"

namespace tinymesh {

class TINYMESH_API Vertex {
public:
    // Forward declaration
    class VertexIterator;
    class ConstVertexIterator;
    class HalfedgeIterator;
    class ConstHalfedgeIterator;
    class FaceIterator;
    class ConstFaceIterator;

public:
    Vertex();
    Vertex(const Vertex &v) = default;
    Vertex(Vertex &&p) noexcept = default;
    explicit Vertex(const Vec3 &v);
    virtual ~Vertex() = default;

    Vertex &operator=(const Vertex &p) = default;
    Vertex &operator=(Vertex &&p) noexcept = default;

    bool operator==(const Vertex &other) const;

    //! Number of connected vertices
    int degree() const;

    //! Compute vertex normal with those of surrounding faces
    Vec3 normal() const;

    //! Obtuse-safe Volonoi area calculation [Meyer et al. 2003]
    double volonoiArea() const;

    //! Volonoi area in specified face [Meyer et al. 2003]
    double volonoiArea(const Face *const f) const;

    /**
     * Gaussian curvature for vertex [Meyer et al. 2003]
     */
    double K() const;
    /**
     * Mean curvature for vertex [Meyer et al. 2003] 
     */
    double H() const;

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

    Vec3 pos() const {
        return pos_;
    }
    void setPos(const Vec3 &pos) {
        pos_ = pos;
    }

    int index() const {
        return index_;
    }

    bool isBoundary() const;

    bool isLocked() const {
        return isLocked_;
    }

    void lock() {
        isLocked_ = true;
    }

    void unlock() {
        isLocked_ = false;
    }

private:
    Vec3 pos_;
    Halfedge *halfedge_ = nullptr;
    int index_ = -1;
    bool isLocked_ = false;

    friend class Mesh;
};

/**
 * VertexIterator
 * @detail Traverse neighboring vertices in the clockwise order.
 */
class TINYMESH_API Vertex::VertexIterator {
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

/**
 * ConstVertexIterator
 * @detail Traverse neighboring vertices in the clockwise order. 
 */
class TINYMESH_API Vertex::ConstVertexIterator {
public:
    explicit ConstVertexIterator(Halfedge *he);
    bool operator!=(const ConstVertexIterator &it) const;
    const Vertex &operator*() const;
    const Vertex *ptr() const;
    const Vertex *operator->() const;
    ConstVertexIterator &operator++();
    ConstVertexIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

/**
 * HalfedgeIterator
 * @detail Traverse outward halfedges in the clockwise order.
 */
class TINYMESH_API Vertex::HalfedgeIterator {
public:
    HalfedgeIterator(Halfedge *he);
    bool operator!=(const HalfedgeIterator &it) const;
    Halfedge &operator*();
    Halfedge *ptr() const;
    Halfedge *operator->() const;
    HalfedgeIterator &operator++();
    HalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

/**
 * ConstHalfedgeIterator
 * @detail Traverse outward halfedges in the clockwise order.
 */

class TINYMESH_API Vertex::ConstHalfedgeIterator {
public:
    ConstHalfedgeIterator(Halfedge *he);
    bool operator!=(const ConstHalfedgeIterator &it) const;
    const Halfedge &operator*() const;
    const Halfedge *ptr() const;
    const Halfedge *operator->() const;
    ConstHalfedgeIterator &operator++();
    ConstHalfedgeIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

/**
 * FaceIterator
 * @detail Traverse neighboring faces in clockwise order.
 */
class TINYMESH_API Vertex::FaceIterator {
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

/**
 * ConstFaceIterator
 * @detail Traverse neighboring faces in clockwise order.
 */
class TINYMESH_API Vertex::ConstFaceIterator {
public:
    ConstFaceIterator(Halfedge *he);
    bool operator!=(const ConstFaceIterator &it) const;
    const Face &operator*() const;
    const Face *ptr() const;
    const Face *operator->() const;
    ConstFaceIterator &operator++();
    ConstFaceIterator operator++(int);

private:
    Halfedge *halfedge_, *iter_;
};

}  // namespace tinymesh

#endif  // TINYMESH_VERTEX_H
