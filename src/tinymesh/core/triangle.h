#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_TRIANGLE_H
#define TINYMESH_TRIANGLE_H

#include "api.h"
#include "vec.h"
#include "bounds.h"

namespace tinymesh {

class TINYMESH_API Triangle {
public:
    Triangle() = default;
    Triangle(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2)
        : p0(p0)
        , p1(p1)
        , p2(p2) {
    }

    double distance(const Vec3 &p) const;
    Vec3 closestPoint(const Vec3 &p) const;
    Bounds3 bounds() const;

private:
    Vec3 p0, p1, p2;
};

}  // namespace tinymesh

#endif  // TINYMESH_TRIANGLE_H
