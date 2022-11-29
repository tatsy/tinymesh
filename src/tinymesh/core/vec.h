#pragma once

#include <array>
#include <algorithm>
#include <functional>
#include <type_traits>

template <typename Float, int Dims>
class Vec {
    static_assert(std::is_floating_point<Float>::value, "Vector base type must be floating point number!");
    static_assert(Dims > 1, "Vector dimension must be more than 1!");

public:
    Vec() {
        std::fill(elems.begin(), elems.end(), Float(0.0));
    }

    explicit Vec(Float x) {
        std::fill(elems.begin(), elems.end(), x);
    }

    explicit Vec(Float x, Float y, Float z = 0.0, Float w = 0.0) {
        std::fill(elems.begin(), elems.end(), w);

        elems[0] = x;
        if (Dims >= 1) elems[1] = y;
        if (Dims >= 2) elems[2] = z;
    }

    bool operator==(const Vec &other) const {
        if (elems.size() != other.elems.size()) {
            return false;
        }

        for (int i = 0; i < Dims; i++) {
            if (elems[i] != other.elems[i]) {
                return false;
            }
        }

        return true;
    }

    bool operator!=(const Vec &other) const {
        return !this->operator==(other);
    }

    Vec &operator+=(const Vec &other) {
        if (elems.size() != other.elems.size()) {
            throw std::runtime_error("Dimensions do not match!");
        }

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
        if (elems.size() != other.elems.size()) {
            throw std::runtime_error("Dimensions do not match!");
        }

        for (int i = 0; i < Dims; i++) {
            if (other.elems[i] == 0.0) {
                throw std::runtime_error("Zero division!");
            }
            elems[i] /= other.elems[i];
        }

        return *this;
    }

    Vec &operator/=(double s) {
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

    template <typename Dummy = Float>
    typename std::enable_if<Dims >= 1, Dummy>::type x() const {
        return elems[0];
    }

    template <typename Dummy = Float>
    typename std::enable_if<Dims >= 2, Dummy>::type y() const {
        return elems[1];
    }

    template <typename Dummy = Float>
    typename std::enable_if<Dims >= 3, Dummy>::type z() const {
        return elems[2];
    }

    template <typename Dummy = Float>
    typename std::enable_if<Dims >= 4, Dummy>::type w() const {
        return elems[3];
    }

private:
    std::array<Float, Dims> elems;
};

// Type definition
using Vec2 = Vec<double, 2>;
using Vec3 = Vec<double, 3>;

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
    auto ret = v;
    ret *= s;
    return ret;
}

template <typename Float, int Dims>
Vec<Float, Dims> operator-(Float s, const Vec<Float, Dims> &v) {
    auto ret = v;
    ret *= s;
    return ret;
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
Vec<Float, Dims> cross(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2,
                       typename std::enable_if<Dims == 3>::type * = 0) {
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
template <typename Float, int Dims>
Float dihedral(const Vec<Float, Dims> &v1, const Vec<Float, Dims> &v2, const Vec<Float, Dims> &v3,
               const Vec<Float, Dims> &v1rev, typename std::enable_if<Dims == 3>::type * = 0) {
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

// Hash
namespace std {

//template <>
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
