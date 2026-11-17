#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_VEC_H
#define TINYMESH_VEC_H

#include <cmath>
#include <array>
#include <algorithm>
#include <ostream>
#include <functional>
#include <type_traits>
#include <stdexcept>

#include <Eigen/Core>

template <typename Float, int Dims>
    requires std::is_floating_point<Float>::value && (Dims >= 2)
class Vec {
    using EigenVectorType = Eigen::Matrix<Float, Dims, 1>;

public:
    Vec() = default;

    explicit Vec(Float x) {
        std::fill(elems.begin(), elems.end(), x);
    }

    explicit Vec(Float x, Float y, Float z = 0.0, Float w = 0.0) {
        std::fill(elems.begin(), elems.end(), w);

        elems[0] = x;
        elems[1] = y;
        if constexpr (Dims >= 3) elems[2] = z;
    }

    explicit Vec(const Vec &other) = default;
    Vec &operator=(const Vec &other) = default;

    Vec(const EigenVectorType &v) {
        for (int d = 0; d < Dims; d++) {
            elems[d] = v(d);
        }
    }

    Vec &operator=(const EigenVectorType &v) {
        for (int d = 0; d < Dims; d++) {
            elems[d] = v(d);
        }
        return *this;
    }

    operator EigenVectorType() const {
        EigenVectorType v;
        for (int d = 0; d < Dims; d++) {
            v(d) = elems[d];
        }
        return v;
    }

    bool operator==(const Vec &other) const = default;

    bool operator!=(const Vec &other) const {
        return !this->operator==(other);
    }

    Vec &operator+=(const Vec &other) {
        for (size_t i = 0; i < Dims; i++) {
            elems[i] += other.elems[i];
        }

        return *this;
    }

    Vec operator-() const {
        Vec<Float, Dims> ret;
        for (size_t i = 0; i < elems.size(); i++) {
            ret.elems[i] = -elems[i];
        }
        return ret;
    }

    Vec &operator-=(const Vec &other) {
        return this->operator+=(-other);
    }

    Vec &operator*=(const Vec &other) {
        if (elems.size() != other.elems.size()) {
            throw std::runtime_error("Dimensions do not match!");
        }

        for (int i = 0; i < Dims; i++) {
            elems[i] *= other.elems[i];
        }

        return *this;
    }

    Vec &operator*=(Float s) {
        for (int i = 0; i < Dims; i++) {
            elems[i] *= s;
        }

        return *this;
    }

    Vec &operator/=(const Vec &other) {
        for (int i = 0; i < Dims; i++) {
            if (other.elems[i] == 0.0) {
                throw std::runtime_error("Zero division!");
            }
            elems[i] /= other.elems[i];
        }

        return *this;
    }

    Vec &operator/=(Float s) {
        if (s == 0.0) {
            throw std::runtime_error("Zero division");
        }

        for (int i = 0; i < Dims; i++) {
            elems[i] /= s;
        }

        return *this;
    }

    Float &operator[](int i) {
        return elems[i];
    }

    Float operator[](int i) const {
        return elems[i];
    }

    Float x() const noexcept {
        return elems[0];
    }

    Float y() const noexcept {
        return elems[1];
    }

    Float z() const noexcept
        requires(Dims >= 3)
    {
        return elems[2];
    }

    Float w() const noexcept
        requires(Dims >= 4)
    {
        return elems[3];
    }

private:
    std::array<Float, Dims> elems{};
};

// Type definition
using Vec2 = Vec<double, 2>;
using Vec3 = Vec<double, 3>;

// standard output
template <typename Float, int Dims>
std::ostream &operator<<(std::ostream &os, const Vec<Float, Dims> &v) {
    os << "[ ";
    for (int d = 0; d < Dims; d++) {
        os << v[d];
        if (d != Dims - 1) os << ", ";
    }
    os << " ]";
    return os;
}

// Basic arithmetics
template <typename Float, int Dims>
Vec<Float, Dims> operator+(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    auto ret = v1;
    ret += v2;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator-(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    auto ret = v1;
    ret -= v2;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator-(const Vec<Float, Dims> &v, Float s) {
    return v - Vec<Float, Dims>(s);
}

template <typename Float, int Dims>
Vec<Float, Dims> operator-(Float s, const Vec<Float, Dims> &v) {
    return Vec<Float, Dims>(s) - v;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator*(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    auto ret = v1;
    ret *= v2;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator*(const Vec<Float, Dims> &v, Float s) {
    auto ret = v;
    ret *= s;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator*(Float s, const Vec<Float, Dims> &v) {
    auto ret = v;
    ret *= s;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator/(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    auto ret = v1;
    ret /= v2;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator/(const Vec<Float, Dims> &v, Float s) {
    auto ret = v;
    ret /= s;
    return ret;
}

template <typename Float, int Rows, int Cols>
Vec<Float, Rows> operator*(const Eigen::Matrix<Float, Rows, Cols> &m, const Vec<Float, Cols> &v) {
    return Vec<Float, Rows>(m * Eigen::Matrix<Float, Cols, 1>(v));
}

// GLSL like vector arithmetics
template <typename Float, int Dims>
Float dot(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    Float ret = 0.0;
    for (int i = 0; i < Dims; i++) {
        ret += v1[i] * v2[i];
    }
    return ret;
}

template <typename Float, int Dims>
    requires(Dims == 3)
Vec<Float, Dims> cross(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2) {
    const Float x = v1[1] * v2[2] - v1[2] * v2[1];
    const Float y = v1[2] * v2[0] - v1[0] * v2[2];
    const Float z = v1[0] * v2[1] - v1[1] * v2[0];
    return Vec<Float, Dims>{ x, y, z };
}

template <typename Float, int Dims>
Float length(const Vec<Float, Dims> &v) {
    return std::sqrt(dot(v, v));
}

template <typename Float, int Dims>
Vec<Float, Dims> normalize(const Vec<Float, Dims> &v) {
    return v / length(v);
}

// Arithmetics for 3D vectors

//! Angle a-b-c in radian (smaller side)
template <typename Float, int Dims>
    requires(Dims == 3)
double angle(const Vec<Float, Dims> &a, const Vec<Float, Dims> &b, const Vec<Float, Dims> &c) {
    const auto e0 = a - b;
    const auto e1 = c - b;
    return std::atan2(length(cross(e0, e1)), dot(e0, e1));
}

//! check if the triangle is obtuse
template <typename Float, int Dims>
    requires(Dims == 3)
bool obtuse(const Vec<Float, Dims> &a, const Vec<Float, Dims> &b, const Vec<Float, Dims> &c) {
    return dot(b - a, c - a) < 0.0 || dot(a - b, c - b) < 0.0 || dot(a - c, b - c) < 0.0;
}

//! cotangent for angle a-b-c
template <typename Float, int Dims>
    requires(Dims == 3)
Float cot(const Vec<Float, Dims> &a, const Vec<Float, Dims> &b, const Vec<Float, Dims> &c) {
    const auto e0 = a - b;
    const auto e1 = c - b;
    return dot(e0, e1) / (length(cross(e0, e1)) + (Float)1.0e-20);
}

template <typename Float, int Dims>
    requires(Dims == 3)
Float dihedral(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2, const Vec<Float, Dims> &v3,
               const Vec<Float, Dims> &v1rev) {
    using V = Vec<Float, Dims>;
    const V e12 = v2 - v1;
    const V e13 = v3 - v1;
    const V n1 = cross(e12, e13);
    const Float l1 = length(n1);
    const V e42 = v2 - v1rev;
    const V e43 = v3 - v1rev;
    const V n4 = cross(e42, e43);
    const Float l4 = length(n4);
    if (l1 == 0.0 || l4 == 0.0) {
        return 0.0;
    }

    const double cosTheta = dot(n1, n4) / (l1 * l4);
    return std::acos(std::max(-1.0, std::min(cosTheta, 1.0)));
}

namespace std {

// Hash
template <typename Float, int Dims>
struct hash<Vec<Float, Dims>> {
    std::size_t operator()(const Vec<Float, Dims> &v) const {
        std::size_t h = 0;
        for (int i = 0; i < Dims; i++) {
            h = std::hash<Float>()(v[i]) ^ (h << 1);
        }
        return h;
    }
};

}  // namespace std

template <typename Float, int Dims>
Vec<Float, Dims> abs(const Vec<Float, Dims> &v) {
    Vec<Float, Dims> ret;
    for (int d = 0; d < Dims; d++) {
        ret[d] = std::abs(v[d]);
    }
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> min(const Vec<Float, Dims> &v0, const Vec<Float, Dims> &v1) {
    Vec<Float, Dims> v;
    for (int d = 0; d < Dims; d++) {
        v[d] = std::min(v0[d], v1[d]);
    }
    return v;
}

template <typename Float, int Dims>
Vec<Float, Dims> max(const Vec<Float, Dims> &v0, const Vec<Float, Dims> &v1) {
    Vec<Float, Dims> v;
    for (int d = 0; d < Dims; d++) {
        v[d] = std::max(v0[d], v1[d]);
    }
    return v;
}

#endif  // TINYMESH_VEC_H
