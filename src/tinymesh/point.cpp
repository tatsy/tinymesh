#define TINYMESH_API_EXPORT
#include "point.h"

namespace tinymesh {

Point::Point() {}
Point::Point(const Vector &v) : m_pos{ v } {}

}  // naemspace tinymesh
