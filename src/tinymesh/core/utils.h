#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_UTILS_H
#define TINYMESH_UTILS_H

#include <map>
#include <tuple>
#include <functional>

#include "vec.h"
#include "eigen.h"
#include "constants.h"

namespace tinymesh {

//! Pair of indices (i.e., unsigned int)
using IndexPair = std::tuple<uint32_t, uint32_t>;

struct IndexPairHash {
    std::size_t operator()(const IndexPair &key) const noexcept {
        const auto h1 = std::hash<std::uint32_t>{}(std::get<0>(key));
        const auto h2 = std::hash<std::uint32_t>{}(std::get<1>(key));

        return h1 ^ (h2 + std::size_t{ 0x9e3779b9u } + (h1 << 6) + (h1 >> 2));
    }
};

//! Local coordinate frame
using LocalFrame = std::tuple<Vec3, Vec3, Vec3>;

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

//! Rotate vector v0 to vector v1
inline EigenMatrix3 rotationFromTwoVectors(const Vec3 &v0, const Vec3 &v1) {
    const Vec3 outer = cross(v0, v1);
    const double inner = dot(v0, v1);
    const double outerNorm = length(outer);
    const double scale = length(v0) * length(v1);
    if (outerNorm < 1.0e-12 * scale) {
        if (inner > 0.0) {
            // No rotation needed
            return EigenMatrix3::Identity();
        } else {
            // 180 degree rotation
            Vec3 axis;
            if (std::abs(v0.x()) < std::abs(v0.y())) {
                axis = Vec3{ 0.0, -v0.z(), v0.y() };
            } else {
                axis = Vec3{ -v0.z(), 0.0, v0.x() };
            }
            axis = normalize(axis);
            return rotationAxisAngle(kPi, axis);
        }
    }

    const double theta = std::atan2(outerNorm, inner);
    const Vec3 axis = normalize(outer);
    return rotationAxisAngle(theta, axis);
}

}  // namespace tinymesh

#endif  // TINYMESH_UTILS_H
