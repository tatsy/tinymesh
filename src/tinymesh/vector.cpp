#define TINYMESH_API_EXPORT
#include "vector.h"

namespace tinymesh {
    
Vector::Vector() {}
Vector::Vector(double x) : x{ x }, y{ x }, z{ x } {}
Vector::Vector(double x, double y, double z) : x{ x }, y{ y }, z{ z } {}
Vector::Vector(const Vector &v) : x{ v.x }, y{ v.y }, z{ v.z } {}

Vector &Vector::operator=(const Vector &v) {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
    return *this;
}

bool Vector::operator==(const Vector &v) const {
    return x == v.x && y == v.y && z == v.z;
}

bool Vector::operator!=(const Vector &v) const {
    return !this->operator==(v);
}

Vector &Vector::operator+=(const Vector &v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

Vector Vector::operator-() const {
    return Vector(-x, -y, -z);
}

Vector &Vector::operator-=(const Vector &v) {
    this->operator+=(-v);
    return *this;
}

Vector &Vector::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

Vector &Vector::operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

double Vector::dot(const Vector &v) const {
    return x * v.x + y * v.y + z * v.z;
}

Vector Vector::cross(const Vector &v) const {
    double xx = y * v.z - z * v.y;
    double yy = z * v.x - x * v.z;
    double zz = x * v.y - y * v.x;
    return Vector(xx, yy, zz);
}

double Vector::length() const {
    return std::sqrt(this->dot(*this));
}

void Vector::normalize() {
    this->operator/=(this->length());
}

}  // namespace tinymesh

tinymesh::Vector operator+(const tinymesh::Vector &v1, const tinymesh::Vector &v2) {
    using tinymesh::Vector;
    Vector v = v1;
    v += v2;
    return v;
}

tinymesh::Vector operator-(const tinymesh::Vector &v1, const tinymesh::Vector &v2) {
    using tinymesh::Vector;
    Vector v = v1;
    v -= v2;
    return v;
}

tinymesh::Vector operator*(const tinymesh::Vector &v, double s) {
    using tinymesh::Vector;
    Vector u = v;
    u *= s;
    return u;
}

tinymesh::Vector operator*(double s, const tinymesh::Vector &v) {
    using tinymesh::Vector;
    Vector u = v;
    u *= s;
    return u;
}

tinymesh::Vector operator/(const tinymesh::Vector &v, double s) {
    using tinymesh::Vector;
    Vector u = v;
    u /= s;
    return u;
}

