#define TINYMESH_API_EXPORT
#include "halfedge.h"

#include "core/vertex.h"
#include "core/face.h"

namespace tinymesh {

Halfedge::Halfedge() {
}

bool Halfedge::operator==(const Halfedge &other) const {
    bool ret = true;
    ret &= (src_ == other.src_);
    ret &= (next_ == other.next_);
    ret &= (rev_ == other.rev_);
    ret &= (edge_ == other.edge_);
    ret &= (face_ == other.face_);
    ret &= (index_ == other.index_);
    return ret;
}

double Halfedge::length() const {
    return ::length(src()->pos() - dst()->pos());
}

double Halfedge::cotWeight() const {
    Vertex *v0 = this->src();
    Vertex *v1 = this->dst();
    Vertex *v2 = this->next()->dst();
    Vertex *v3 = this->rev()->next()->dst();
    const Vec3 &p0 = v0->pos();
    const Vec3 &p1 = v1->pos();
    const Vec3 &p2 = v2->pos();
    const Vec3 &p3 = v3->pos();
    const double sin_a = ::length(cross(p0 - p2, p1 - p2));
    const double cos_a = dot(p0 - p2, p1 - p2);
    const double cot_a = cos_a / std::max(sin_a, 1.0e-6);
    const double sin_b = ::length(cross(p0 - p3, p1 - p3));
    const double cos_b = dot(p0 - p3, p1 - p3);
    const double cot_b = cos_b / std::max(sin_b, 1.0e-6);
    return 0.5 * (cot_a + cot_b);
}

bool Halfedge::isLocked() const {
    return this->src()->isLocked() && this->dst()->isLocked();
}

bool Halfedge::isBoundary() const {
    return this->src()->isBoundary() && this->dst()->isBoundary();
}

}  // namespace tinymesh
