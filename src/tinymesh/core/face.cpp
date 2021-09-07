#define TINYMESH_API_EXPORT
#include "face.h"

#include "vertex.h"
#include "halfedge.h"

namespace tinymesh {

Face::Face() {
}

bool Face::operator==(const Face &other) const {
    bool ret = true;
    ret &= (halfedge_ == other.halfedge_);
    ret &= (index_ == other.index_);
    return ret;
}

bool Face::isBoundary() {
    for (auto vit = this->v_begin(); vit != this->v_end(); ++vit) {
        if (vit->isBoundary()) {
            return true;
        }
    }
    return false;
}

bool Face::isStatic() {
    for (auto vit = this->v_begin(); vit != this->v_end(); ++vit) {
        if (vit->isStatic()) {
            return true;
        }
    }
    return false;
}

Face::VertexIterator Face::v_begin() {
    return Face::VertexIterator(halfedge_);
}

Face::VertexIterator Face::v_end() {
    return Face::VertexIterator(nullptr);
}

Face::FaceIterator Face::f_begin() {
    return Face::FaceIterator(halfedge_);
}

Face::FaceIterator Face::f_end() {
    return Face::FaceIterator(nullptr);
}

// ----------
// VertexIterator
// ----------

Face::VertexIterator::VertexIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::VertexIterator::operator!=(const Face::VertexIterator &it) const {
    return iter_ != it.iter_;
}

Vertex &Face::VertexIterator::operator*() {
    return *iter_->src();
}

Vertex *Face::VertexIterator::ptr() const {
    return iter_->src();
}

Vertex *Face::VertexIterator::operator->() const {
    return iter_->src();
}

Face::VertexIterator &Face::VertexIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::VertexIterator Face::VertexIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::VertexIterator(tmp);
}

// ----------
// FaceIterator
// ----------

Face::FaceIterator::FaceIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::FaceIterator::operator!=(const Face::FaceIterator &it) const {
    return iter_ != it.iter_;
}

Face &Face::FaceIterator::operator*() {
    return *iter_->rev()->face();
}

Face *Face::FaceIterator::ptr() const {
    return iter_->rev()->face();
}

Face *Face::FaceIterator::operator->() const {
    return iter_->rev()->face();
}

Face::FaceIterator &Face::FaceIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::FaceIterator Face::FaceIterator::operator++(int) {
    Halfedge* tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::FaceIterator(tmp);
}

}  // namespace tinymesh
