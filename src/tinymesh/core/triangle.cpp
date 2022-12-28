#define TINYMESH_API_EXPORT
#include "triangle.h"

namespace tinymesh {

double Triangle::distance(const Vec3 &p) const {
    return length(p - closestPoint(p));
}

Vec3 Triangle::closestPoint(const Vec3 &p) const {
    const Vec3 ab = p1 - p0;
    const Vec3 ac = p2 - p0;
    const Vec3 ap = p - p0;

    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return p0;  //#1

    const Vec3 bp = p - p2;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return p1;  //#2

    const Vec3 cp = p - p2;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return p2;  //#3

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return p0 + v * ab;  //#4
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double v = d2 / (d2 - d6);
        return p0 + v * ac;  //#5
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return p1 + v * (p2 - p1);  //#6
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return p0 + v * ab + w * ac;  //#0
}

Bounds3 Triangle::bounds() const {
    Bounds3 ret;
    ret.merge(p0);
    ret.merge(p1);
    ret.merge(p2);
    return ret;
}

}  // namespace tinymesh
