#pragma once

#include <functional>

#include "core/common.h"

class TINYMESH_API Vec {
public:
    Vec();
    explicit Vec(double x);
    Vec(double x, double y, double z);

    bool operator==(const Vec &other) const;
    bool operator!=(const Vec &other) const;

    Vec &operator+=(const Vec &other);
    Vec operator-() const;
    Vec &operator-=(const Vec &other);
    Vec &operator*=(const Vec &other);
    Vec &operator*=(double s);
    Vec &operator/=(const Vec &other);
    Vec &operator/=(double s);

    double operator[](int i) const;

    double x, y, z;
};

// Basic arithmetics
TINYMESH_API Vec operator-(const Vec &v1, const Vec &v2);
TINYMESH_API Vec operator*(const Vec &v1, const Vec &v2);
TINYMESH_API Vec operator+(const Vec &v1, const Vec &v2);
TINYMESH_API Vec operator*(const Vec &v1, double s);
TINYMESH_API Vec operator*(double s, const Vec &v2);
TINYMESH_API Vec operator/(const Vec &v1, const Vec &v2);
TINYMESH_API Vec operator/(const Vec &v1, double s);

// GLSL like vector arithmetics
TINYMESH_API double dot(const Vec &v1, const Vec &v2);
TINYMESH_API Vec cross(const Vec &v1, const Vec &v2);
TINYMESH_API Vec normalize(const Vec &v);
TINYMESH_API double length(const Vec &v);

// Hash
namespace std {

template <>
struct hash<Vec> {
    std::size_t operator()(const Vec& v) const {
        std::size_t h = 0;
        h = std::hash<double>()(v.x) ^ (h << 1);
        h = std::hash<double>()(v.y) ^ (h << 1);
        h = std::hash<double>()(v.z) ^ (h << 1);
        return h;
    }
};

}  // namespace std
