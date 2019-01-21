#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_HALFEDGE_H
#define TINYMESH_HALFEDGE_H

#include <memory>

#include "common.h"

namespace tinymesh {

class Halfedge {
public:
    Halfedge();

private:
    std::weak_ptr<Point> m_start, m_end;
    std::weak_ptr<Halfedge> m_next;
    std::weak_ptr<Face> m_face;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_HALFEDGE_H
