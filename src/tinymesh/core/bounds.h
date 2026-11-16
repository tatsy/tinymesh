#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_BOUNDS_H
#define TINYMESH_BOUNDS_H

#include "debug.h"
#include "vec.h"

namespace tinymesh {

template <typename Float, int Dims>
    requires std::is_floating_point<Float>::value && (Dims >= 2)
class Bounds {
    typedef Vec<Float, Dims> VecD;

public:
    Bounds() = default;

    Bounds(const VecD &posMin, const VecD &posMax)
        : posMin_(posMin)
        , posMax_(posMax) {
    }

    void merge(const VecD &p) {
        posMin_ = ::min(posMin_, p);
        posMax_ = ::max(posMax_, p);
    }

    double area() const
        requires(Dims == 3)
    {
        const VecD d = ::abs(posMax_ - posMin_);
        return (2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]));
    }

    bool inside(const VecD &p) const {
        for (int d = 0; d < Dims; d++) {
            if (p[d] < posMin_[d] || posMax_[d] < p[d]) return false;
        }
        return true;
    }

    int maxExtent() const {
        const VecD v = ::abs(posMax_ - posMin_);
        int dim = -1;
        Float vmax = -1.0;
        for (int d = 0; d < Dims; d++) {
            if (vmax < v[d]) {
                vmax = v[d];
                dim = d;
            }
        }

        Assertion(dim >= 0, "Something is wrong!");
        return dim;
    }

    double distance(const VecD &p) {
        const Vec3 vmin = ::abs(posMin_ - p);
        const Vec3 vmax = ::abs(posMax_ - p);
        double dist = 0.0;
        for (int d = 0; d < Dims; d++) {
            if (p[d] < posMin_[d] || posMax_[d] < p[d]) {
                const double gap = std::min(vmin[d], vmax[d]);
                dist += gap * gap;
            }
        }
        return std::sqrt(dist);
    }

    static Bounds<Float, Dims> merge(const Bounds<Float, Dims> &b0, const Bounds<Float, Dims> &b1) {
        Bounds<Float, Dims> ret;
        ret.posMin_ = ::min(b0.posMin_, b1.posMin_);
        ret.posMax_ = ::max(b0.posMax_, b1.posMax_);
        return ret;
    }

    VecD posMin() const {
        return posMin_;
    }
    VecD posMax() const {
        return posMax_;
    }

private:
    VecD posMin_ = VecD(std::numeric_limits<Float>::max());
    VecD posMax_ = VecD(std::numeric_limits<Float>::lowest());
};

using Bounds2 = Bounds<double, 2>;
using Bounds3 = Bounds<double, 3>;

}  // namespace tinymesh

#endif  // TINYMESH_BOUNDS_H
