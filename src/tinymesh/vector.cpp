#define TINYMESH_API_EXPORT
#include "vector.h"

namespace tinymesh {
    
Vector::Vector() {}
Vector::Vector(double x) : x{ x }, y{ x }, z{ x } {}
Vector::Vector(double x, double y, double z) : x{ x }, y{ y }, z{ z } {}

}  // namespace tinymesh
