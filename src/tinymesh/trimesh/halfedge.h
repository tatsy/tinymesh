#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_HALFEDGE_H
#define TINYMESH_HALFEDGE_H

#include <memory>

#include "core/common.h"

namespace tinymesh {

class TINYMESH_API Halfedge {
public:
    Halfedge();
    Halfedge(const Halfedge &he);
    Halfedge(Halfedge &&he) noexcept;
    virtual ~Halfedge() = default;

    double length() const;
    bool isBorder() { return rev_ == nullptr; }

    Halfedge &operator=(const Halfedge &he) = default;
    Halfedge &operator=(Halfedge &&he) noexcept = default;

    Vertex *src() const { return src_; }
    Vertex *dst() const { return next_->src_; }
    Halfedge *next() const { return next_; }
    Halfedge *prev() const { return next_->next_; }
    Halfedge *rev() const { return rev_; }
    Face *face() const { return face_; }
    int index() const { return index_; }

private:
    Vertex *src_ = nullptr;
    Halfedge *next_ = nullptr;
    Halfedge *rev_ = nullptr;
    Face *face_ = nullptr;
    int index_ = -1;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_HALFEDGE_H
