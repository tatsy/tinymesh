#define TINYMESH_API_EXPORT
#include "halfedge.h"

#include "vertex.h"
#include "face.h"

namespace tinymesh {

Halfedge::Halfedge() {}

Halfedge::Halfedge(const tinymesh::Halfedge &he)
    : Halfedge{} {
    this->operator=(he);
}

Halfedge::Halfedge(tinymesh::Halfedge &&he) noexcept
    : Halfedge{} {
    this->operator=(std::move(he));
}

Halfedge &Halfedge::operator=(const Halfedge &he) {
    src_ = he.src_;
    next_ = he.next_;
    rev_ = he.rev_;
    face_ = he.face_;
    index_ = he.index_;
    return *this;
}

Halfedge &Halfedge::operator=(Halfedge &&he) noexcept {
    src_ = he.src_;
    next_ = he.next_;
    rev_ = he.rev_;
    face_ = he.face_;
    index_ = he.index_;

    he.src_ = nullptr;
    he.next_ = nullptr;
    he.rev_ = nullptr;
    he.face_ = nullptr;
    he.index_ = -1;
    return *this;
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

bool Halfedge::isBoundary() const {
    return face_->isBoundary();
}

}  // namespace tinymesh
