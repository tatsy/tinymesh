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

private:
    std::weak_ptr<Halfedge> m_he;

    friend class Mesh;
};

}  // namespace tinymesh

#endif  // TINYMESH_FACE_H
