#ifdef _MSC_VER
#pragma once
#endif 

#ifndef TINYMESH_FACE_H
#define TINYMESH_FACE_H

#include <memory>

#include "common.h"

namespace tinymesh {

class TINYMESH_EXPORTS Face {
public:
    Face();
    Face(const Face &face) = default;
    Face(Face &&face) noexcept = default;
    virtual ~Face() = default;

    Face &operator=(const Face &face) = default;
    Face &operator=(Face &&face) noexcept = default;

private:
    Halfedge *m_he;
    int m_index;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_FACE_H
