#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_VECTOR_H
#define TINYMESH_VECTOR_H

#include "common.h"

namespace tinymesh {
    
class TINYMESH_EXPORTS Vector {
public:
    Vector();
    explicit Vector(double x);
    Vector(double x, double y, double z);

    double x, y, z;
};

}  // namespace tinymesh

#endif  // TINYMESH_VECTOR_H
