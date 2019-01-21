#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_POINT_H
#define TINYMESH_POINT_H

#include <memory>

#include "common.h"
#include "vector.h"

namespace tinymesh {

class TINYMESH_EXPORTS Point {
public:
    Point();
    Point(const Vector &v);
    inline Vector pos() const { return m_pos; }

private:
    Vector m_pos;
    std::weak_ptr<Halfedge> m_he;
    std::weak_ptr<Face> m_face;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_POINT_H
