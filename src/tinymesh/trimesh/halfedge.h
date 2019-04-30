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
    bool isBorder() { return m_rev == nullptr; }

    Halfedge &operator=(const Halfedge &he) = default;
    Halfedge &operator=(Halfedge &&he) noexcept = default;

    Vertex *src() const { return m_src; }
    Vertex *dst() const { return m_next->m_src; }
    Halfedge *next() const { return m_next; }
    Halfedge *prev() const { return m_next->m_next; }
    Halfedge *rev() const { return m_rev; }
    Face *face() const { return m_face; }
    int index() const { return m_index; }

private:
    Vertex *m_src = nullptr;
    Halfedge *m_next = nullptr;
    Halfedge *m_rev = nullptr;
    Face *m_face = nullptr;
    int m_index;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_HALFEDGE_H
