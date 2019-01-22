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

    bool operator==(const Vector &v) const;
    bool operator!=(const Vector &v) const;

    Vector &operator+=(const Vector &v);
    Vector operator-() const;
    Vector &operator-=(const Vector &v);
    Vector &operator*=(double s);
    Vector &operator/=(double s);

    double dot(const Vector &v) const;
    Vector cross(const Vector &v) const;
    double length() const;

    void normalize();

    double x, y, z;
};

}  // namespace tinymesh

TINYMESH_EXPORTS tinymesh::Vector operator+(const tinymesh::Vector &v1, const tinymesh::Vector &v2);
TINYMESH_EXPORTS tinymesh::Vector operator-(const tinymesh::Vector &v1, const tinymesh::Vector &v2);
TINYMESH_EXPORTS tinymesh::Vector operator*(const tinymesh::Vector &v, double s);
TINYMESH_EXPORTS tinymesh::Vector operator*(double s, const tinymesh::Vector &v);

#endif  // TINYMESH_VECTOR_H
