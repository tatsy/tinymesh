#include "halfedge.h"

#include <algorithm>

#include "core/vertex.h"
#include "core/face.h"

namespace tinymesh {

namespace {

double cotangent(const Halfedge *he) {
    if (he->face()->isHole()) {
        return 0.0;
    }

    const Vec3 p0 = he->src()->pos();
    const Vec3 p1 = he->dst()->pos();
    const Vec3 p2 = he->next()->dst()->pos();
    const Vec3 e0 = p0 - p2;
    const Vec3 e1 = p1 - p2;

    const double denominator = ::length(cross(e0, e1));
    return dot(e0, e1) / std::max(denominator, 1.0e-6);
}

}  // anonymous namespace

Halfedge::Halfedge() {
}

bool Halfedge::operator==(const Halfedge &other) const {
    bool ret = true;
    ret &= (src_ == other.src_);
    ret &= (next_ == other.next_);
    ret &= (rev_ == other.rev_);
    ret &= (face_ == other.face_);
    ret &= (index_ == other.index_);
    return ret;
}

double Halfedge::length() const {
    return ::length(src()->pos() - dst()->pos());
}

double Halfedge::cotWeight() const {
    return 0.5 * (cotangent(this) + cotangent(rev()));
}

bool Halfedge::isLocked() const {
    return this->src()->isLocked() && this->dst()->isLocked();
}

}  // namespace tinymesh
