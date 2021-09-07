#define TINYMESH_API_EXPORT
#include "halfedge.h"

#include "vertex.h"
#include "face.h"

namespace tinymesh {

Halfedge::Halfedge() {}

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

bool Halfedge::isStatic() const {
    return this->src()->isStatic() && this->dst()->isStatic();
}

bool Halfedge::isBoundary() const {
    return this->src()->isBoundary() && this->dst()->isBoundary();
}

}  // namespace tinymesh
