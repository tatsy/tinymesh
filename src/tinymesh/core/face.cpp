#define TINYMESH_API_EXPORT
#include "face.h"

#include <vector>

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

Vec3 Face::normal() {
    std::vector<Vec3> vs;
    for (auto it = v_begin(); it != v_end(); ++it) {
        vs.push_back(it->pos());
    }

    const int N = static_cast<int>(vs.size());
    Vec3 norm(0.0);
    for (int i = 0; i < N; i++) {
        const int prev = (i - 1 + N) % N;
        const int next = (i + 1) % N;
        const Vec3 &p0 = vs[prev];
        const Vec3 &p1 = vs[i];
        const Vec3 &p2 = vs[next];
        norm += cross(p2 - p1, p0 - p1);
    }
    return normalize(norm);
}

double Face::area() {
    std::vector<Vec3> vs;
    for (auto it = v_begin(); it != v_end(); ++it) {
        vs.push_back(it->pos());
    }

    const int N = static_cast<int>(vs.size());
    double area = 0.0;
    for (int i = 1; i < N - 1; i++) {
        const Vec3 &p0 = vs[i];
        const Vec3 &p1 = vs[i - 1];
        const Vec3 &p2 = vs[i + 1];
        area += 0.5 * length(cross(p1 - p0, p2 - p0));
    }
    return area;
}

bool Face::isHole() {
    // Face is hole if all the vertices are at the boundary.
    for (auto it = this->he_begin(); it != this->he_end(); ++it) {
        if (!it->isBoundary()) {
            return false;
        }
    }
    return true;
}

bool Face::isBoundary() {
    for (auto it = this->he_begin(); it != this->he_end(); ++it) {
        if (it->isBoundary()) {
            return true;
        }
    }
    return false;
}

bool Face::isLocked() {
    for (auto it = this->v_begin(); it != this->v_end(); ++it) {
        if (it->isLocked()) {
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

Face::HalfedgeIterator Face::he_begin() {
    return Face::HalfedgeIterator(halfedge_);
}

Face::HalfedgeIterator Face::he_end() {
    return Face::HalfedgeIterator(nullptr);
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
// HalfedgeIterator
// ----------

Face::HalfedgeIterator::HalfedgeIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::HalfedgeIterator::operator!=(const Face::HalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

Halfedge &Face::HalfedgeIterator::operator*() {
    return *iter_;
}

Halfedge *Face::HalfedgeIterator::ptr() const {
    return iter_;
}

Halfedge *Face::HalfedgeIterator::operator->() const {
    return iter_;
}

Face::HalfedgeIterator &Face::HalfedgeIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::HalfedgeIterator Face::HalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::HalfedgeIterator(tmp);
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
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::FaceIterator(tmp);
}

}  // namespace tinymesh
