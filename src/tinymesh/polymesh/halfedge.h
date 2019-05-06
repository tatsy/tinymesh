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
    Halfedge(const Halfedge &he) = default;
    Halfedge(Halfedge &&he) noexcept = default;
    virtual ~Halfedge() = default;

    Halfedge &operator=(const Halfedge &he) = default;
    Halfedge &operator=(Halfedge &&he) noexcept = default;

    bool operator==(const Halfedge &other) const;

    double length() const;

    Vertex *src() const { return src_; }
    Vertex *dst() const { return next_->src_; }
    Halfedge *next() const { return next_; }
    Halfedge *prev() const {
        Halfedge *iter = next_;
        while (iter->next_ != this) {
            iter = iter->next_;
        }
        return iter;
    }
    Halfedge *rev() const { return rev_; }
    Face *face() const { return face_; }
    int index() const { return index_; }
    bool isBoundary() const;

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
