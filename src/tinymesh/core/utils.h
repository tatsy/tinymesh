#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_UTILS_H
#define TINYMESH_UTILS_H

#include <map>
#include <functional>

#include "vec.h"
#include "eigen.h"

//! Pair of indices (i.e., unsigned int)
using IndexPair = std::pair<uint32_t, uint32_t>;

namespace std {

//! Hash type for IndexPair
template <>
struct hash<IndexPair> {
    std::size_t operator()(const IndexPair &k) const {
        return std::get<0>(k) ^ std::get<1>(k);
    }
};

}  // namespace std

//! Rodrigues rotation formula
inline EigenMatrix3 matrixCrossProd(const Vec3 &w) {
    EigenMatrix3 K;
    K << 0.0, -w.z(), w.y(),  // 1st row
        w.z(), 0.0, -w.x(),   // 2nd row
        -w.y(), w.x(), 0.0;   // 3rd row
    return K;
}

inline EigenMatrix3 rotationAxisAngle(double theta, const Vec3 &axis) {
    EigenMatrix3 I = EigenMatrix3::Identity();
    EigenMatrix3 K = matrixCrossProd(axis);
    return I + std::sin(theta) * K + (1.0 - std::cos(theta)) * K * K;
}

#endif  // TINYMESH_UTILS_H
