#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_EDGE_H
#define TINYMESH_EDGE_H

#include "core/common.h"

namespace tinymesh {

class TINYMESH_API Edge {
public:
    Edge() = default;
    Edge(const Edge &) = default;
    Edge(Edge &&) noexcept = default;
    virtual ~Edge() = default;

    Edge &operator=(const Edge &) = default;
    Edge &operator=(Edge &&) noexcept = default;

    bool operator==(const Edge &other) const;

    Halfedge *halfedge() const { return halfedge_; }

    int index() const { return index_; }
    bool isBoundary() const;
    double length() const;

private:
    Halfedge *halfedge_ = nullptr;
    int index_ = -1;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_EDGE_H
