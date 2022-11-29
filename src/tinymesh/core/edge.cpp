#define TINYMESH_API_EXPORT
#include "edge.h"

#include "core/vertex.h"
#include "core/halfedge.h"

namespace tinymesh {

bool Edge::isBoundary() const {
    return halfedge_->isBoundary() || halfedge_->rev()->isBoundary();
}

bool Edge::isLocked() const {
    return halfedge_->isLocked();
}

bool Edge::operator==(const Edge &other) const {
    bool ret = true;
    ret &= (halfedge_ == other.halfedge_);
    ret &= (index_ == other.index_);
    return ret;
}

double Edge::length() const {
    return halfedge_->length();
}

double Edge::cotWeight() const {
    return halfedge_->cotWeight();
}

}  // namespace tinymesh
