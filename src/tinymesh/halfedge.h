#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_HALFEDGE_H
#define TINYMESH_HALFEDGE_H

#include <memory>

#include "common.h"

namespace tinymesh {

class TINYMESH_EXPORTS Halfedge {
public:
    Halfedge();
    Halfedge(const Halfedge &he);
    Halfedge(Halfedge &&he) noexcept;
    virtual ~Halfedge() = default;

    Halfedge &operator=(const Halfedge &he) = default;
    Halfedge &operator=(Halfedge &&he) noexcept = default;

    Vertex *v_start() const { return m_start; }
    Vertex *v_end() const { return m_end; }
    Halfedge *next() const { return m_next; }
    Halfedge *prev() const { return m_prev; }
    Halfedge *opposite() const { return m_opposite; }

private:
    Vertex *m_start;
    Vertex *m_end;
    Halfedge *m_next;
    Halfedge *m_prev;
    Halfedge *m_opposite;
    Face *m_face;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_HALFEDGE_H
