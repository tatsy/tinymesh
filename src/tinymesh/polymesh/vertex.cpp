#define TINYMESH_API_EXPORT
#include "vertex.h"

#include <set>

#include "face.h"
#include "halfedge.h"

namespace tinymesh {

// ----------
// Vertex
// ----------

Vertex::Vertex()
    : pos_{}
    , halfedge_{ nullptr }
    , index_{ -1 } {
}

Vertex::Vertex(const Vec3 &pos)
    : pos_{ pos }
    , halfedge_{ nullptr }
    , index_{ -1 } {
}

bool Vertex::operator==(const Vertex &other) const {
    bool ret = true;
    ret &= (pos_ == other.pos_);
    ret &= (halfedge_ == other.halfedge_);
    ret &= (index_ == other.index_);
    return ret;
}

int Vertex::degree() {
    int deg = 0;
    for (auto it = ohe_begin(); it != ohe_end(); ++it) {
        if (!it->face()->isBoundary()) {
            deg++;
        }
    }
    return deg;
}

Vec3 Vertex::normal() {
    Vec3 norm(0.0);
    Vertex *prev = nullptr;
    for (auto it = v_begin(); it != v_end(); ++it) {
        if (prev) {
            const Vec3 e1 = prev->pos() - pos();
            const Vec3 e2 = it->pos() - pos();
            norm += cross(e1, e2);
        }
        prev = it.ptr();
    }
    const Vec3 e1 = prev->pos() - pos();
    const Vec3 e2 = v_begin()->pos() - pos();
    norm += cross(e1, e2);

    return normalize(norm);
}

bool Vertex::isBoundary() {
    for (auto it = ohe_begin(); it != ohe_end(); ++it) {
        if (it->isBoundary()) {
            return true;
        }
    }
    return false;
}

Vertex::VertexIterator Vertex::v_begin() {
    return Vertex::VertexIterator(halfedge_);
}

Vertex::VertexIterator Vertex::v_end() {
    return Vertex::VertexIterator(nullptr);
}

Vertex::InHalfedgeIterator Vertex::ihe_begin() {
    return Vertex::InHalfedgeIterator(halfedge_->rev());
}

Vertex::InHalfedgeIterator Vertex::ihe_end() {
    return Vertex::InHalfedgeIterator(nullptr);
}

Vertex::OutHalfedgeIterator Vertex::ohe_begin() {
    return Vertex::OutHalfedgeIterator(halfedge_);
}

Vertex::OutHalfedgeIterator Vertex::ohe_end() {
    return Vertex::OutHalfedgeIterator(nullptr);
}

Vertex::FaceIterator Vertex::f_begin() {
    return Vertex::FaceIterator(halfedge_);
}

Vertex::FaceIterator Vertex::f_end() {
    return Vertex::FaceIterator(nullptr);
}

// ----------
// VertexIterator
// ----------

Vertex::VertexIterator::VertexIterator(tinymesh::Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::VertexIterator::operator!=(const Vertex::VertexIterator &it) const {
    return iter_ != it.iter_;
}

Vertex &Vertex::VertexIterator::operator*() {
    return *iter_->dst();
}

Vertex* Vertex::VertexIterator::ptr() const {
    return iter_->dst();
}

Vertex *Vertex::VertexIterator::operator->() const {
    return iter_->dst();
}

Vertex::VertexIterator &Vertex::VertexIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::VertexIterator Vertex::VertexIterator::operator++(int) {
    Halfedge *tmp = iter_;
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::VertexIterator(tmp);
}

// ----------
// InHalfedgeIterator
// ----------

Vertex::InHalfedgeIterator::InHalfedgeIterator(tinymesh::Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::InHalfedgeIterator::operator!=(const Vertex::InHalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

Halfedge &Vertex::InHalfedgeIterator::operator*() {
    return *iter_;
}

Halfedge *Vertex::InHalfedgeIterator::ptr() const {
    return iter_;
}

Halfedge *Vertex::InHalfedgeIterator::operator->() const {
    return iter_;
}

Vertex::InHalfedgeIterator &Vertex::InHalfedgeIterator::operator++() {
    iter_ = iter_->next()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::InHalfedgeIterator Vertex::InHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::InHalfedgeIterator(tmp);
}

// ----------
// OutHalfedgeIterator
// ----------

Vertex::OutHalfedgeIterator::OutHalfedgeIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::OutHalfedgeIterator::operator!=(const Vertex::OutHalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

Halfedge &Vertex::OutHalfedgeIterator::operator*() {
    return *iter_;
}

Halfedge *Vertex::OutHalfedgeIterator::ptr() const {
    return iter_;
}

Halfedge *Vertex::OutHalfedgeIterator::operator->() const {
    return iter_;
}

Vertex::OutHalfedgeIterator &Vertex::OutHalfedgeIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::OutHalfedgeIterator Vertex::OutHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::OutHalfedgeIterator(tmp);
}

// ----------
// FaceIterator
// ----------

Vertex::FaceIterator::FaceIterator(tinymesh::Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::FaceIterator::operator!=(const Vertex::FaceIterator &it) const {
    return iter_ != it.iter_;
}

Face &Vertex::FaceIterator::operator*() {
    return *iter_->face();
}

Face *Vertex::FaceIterator::ptr() const {
    return iter_->face();
}

Face *Vertex::FaceIterator::operator->() const {
    return iter_->face();
}

Vertex::FaceIterator &Vertex::FaceIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::FaceIterator Vertex::FaceIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::FaceIterator(tmp);
}

}  // naemspace tinymesh
