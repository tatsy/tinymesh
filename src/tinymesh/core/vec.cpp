#define TINYMESH_API_EXPORT
#include "vec.h"

#include <cmath>

#include "core/common.h"

Vec::Vec() 
    : x{ 0.0 }
    , y{ 0.0 }
    , z{ 0.0 } {
}

Vec::Vec(double x)
    : x{ x }
    , y{ x }
    , z{ x } {
}

Vec::Vec(double x, double y, double z)
    : x{ x }
    , y{ y }
    , z{ z } {
}

bool Vec::operator==(const Vec& other) const {
    return x == other.x && y == other.y && z == other.z;
}

bool Vec::operator!=(const Vec& other) const {
    return x != other.x || y != other.y || z != other.z;
}

Vec& Vec::operator+=(const Vec& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec Vec::operator-() const {
    return Vec(-x, -y, -z);
}

Vec& Vec::operator-=(const Vec& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec& Vec::operator*=(const Vec& other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
}

Vec& Vec::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

Vec& Vec::operator/=(const Vec& other) {
    if (other.x == 0.0 || other.y == 0.0 || other.z == 0.0) {
        FatalError("Zero division!");
    }

    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
}

Vec& Vec::operator/=(double s) {
    if (s == 0.0) {
        FatalError("Zero division!");
    }

    x /= s;
    y /= s;
    z /= s;
    return *this;
}

double Vec::operator[](int i) const {
    FatalError("Vec subscription more than 2 is specified!");
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    return 0.0;
}

Vec operator+(const Vec& v1, const Vec& v2) {
    Vec r = v1;
    r += v2;
    return r;
}

Vec operator-(const Vec& v1, const Vec& v2) {
    Vec r = v1;
    r -= v2;
    return r;
}

Vec operator*(const Vec& v1, const Vec& v2) {
    Vec r = v1;
    r *= v2;
    return r;
}

Vec operator*(const Vec& v1, double s) {
    Vec r = v1;
    r *= s;
    return r;
}

Vec operator*(double s, const Vec& v2) {
    Vec r = v2;
    r *= s;
    return r;
}

Vec operator/(const Vec& v1, const Vec& v2) {
    Vec r = v1;
    r /= v2;
    return r;
}

Vec operator/(const Vec& v1, double s) {
    Vec r = v1;
    r /= s;
    return r;
}

double dot(const Vec& v1, const Vec& v2) {
    return v1.x* v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec cross(const Vec& v1, const Vec& v2) {
    const double x = v1.y * v2.z - v1.z * v2.y;
    const double y = v1.z * v2.x - v1.x * v2.z;
    const double z = v1.x * v2.y - v1.y * v2.x;
    return Vec(x, y, z);
}

Vec normalize(const Vec& v) {
    return v / length(v);
}

double length(const Vec& v) {
    return std::sqrt(dot(v, v));
}
