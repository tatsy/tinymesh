#define TINYMESH_API_EXPORT
#include "face.h"

#include <vector>

#include "debug.h"
#include "vertex.h"
#include "halfedge.h"
#include "triangle.h"

namespace tinymesh {

Face::Face() {
}

bool Face::operator==(const Face &other) const {
    bool ret = true;
    ret &= (halfedge_ == other.halfedge_);
    ret &= (index_ == other.index_);
    return ret;
}

Triangle Face::toTriangle() const {
    std::vector<Vec3> vs;
    for (auto it = v_begin(); it != v_end(); ++it) {
        vs.push_back(it->pos());
    }

    Assertion(vs.size() == 3, "Non-triangle face cannot be converted to Triangle!");

    return { vs[0], vs[1], vs[2] };
}

Vec3 Face::normal() const {
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

double Face::area() const {
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

int Face::numCorners() const {
    int count = 0;
    for (auto it = v_begin(); it != v_end(); ++it) {
        count++;
    }
    return count;
}

bool Face::isHole() const {
    // Face is hole if all the vertices are at the boundary.
    for (auto it = this->he_begin(); it != this->he_end(); ++it) {
        if (!it->isBoundary()) {
            return false;
        }
    }
    return true;
}

bool Face::isBoundary() const {
    for (auto it = this->he_begin(); it != this->he_end(); ++it) {
        if (it->isBoundary()) {
            return true;
        }
    }
    return false;
}

bool Face::isLocked() const {
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

Face::ConstVertexIterator Face::v_begin() const {
    return Face::ConstVertexIterator(halfedge_);
}

Face::ConstVertexIterator Face::v_end() const {
    return Face::ConstVertexIterator(nullptr);
}

Face::HalfedgeIterator Face::he_begin() {
    return Face::HalfedgeIterator(halfedge_);
}

Face::HalfedgeIterator Face::he_end() {
    return Face::HalfedgeIterator(nullptr);
}

Face::ConstHalfedgeIterator Face::he_begin() const {
    return Face::ConstHalfedgeIterator(halfedge_);
}

Face::ConstHalfedgeIterator Face::he_end() const {
    return Face::ConstHalfedgeIterator(nullptr);
}

Face::FaceIterator Face::f_begin() {
    return Face::FaceIterator(halfedge_);
}

Face::FaceIterator Face::f_end() {
    return Face::FaceIterator(nullptr);
}

Face::ConstFaceIterator Face::f_begin() const {
    return Face::ConstFaceIterator(halfedge_);
}

Face::ConstFaceIterator Face::f_end() const {
    return Face::ConstFaceIterator(nullptr);
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
// ConstVertexIterator
// ----------

Face::ConstVertexIterator::ConstVertexIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::ConstVertexIterator::operator!=(const Face::ConstVertexIterator &it) const {
    return iter_ != it.iter_;
}

const Vertex &Face::ConstVertexIterator::operator*() const {
    return *iter_->src();
}

const Vertex *Face::ConstVertexIterator::ptr() const {
    return iter_->src();
}

const Vertex *Face::ConstVertexIterator::operator->() const {
    return iter_->src();
}

Face::ConstVertexIterator &Face::ConstVertexIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::ConstVertexIterator Face::ConstVertexIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::ConstVertexIterator(tmp);
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
// ConstHalfedgeIterator
// ----------

Face::ConstHalfedgeIterator::ConstHalfedgeIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::ConstHalfedgeIterator::operator!=(const Face::ConstHalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

const Halfedge &Face::ConstHalfedgeIterator::operator*() const {
    return *iter_;
}

const Halfedge *Face::ConstHalfedgeIterator::ptr() const {
    return iter_;
}

const Halfedge *Face::ConstHalfedgeIterator::operator->() const {
    return iter_;
}

Face::ConstHalfedgeIterator &Face::ConstHalfedgeIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::ConstHalfedgeIterator Face::ConstHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::ConstHalfedgeIterator(tmp);
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

// ----------
// ConstFaceIterator
// ----------

Face::ConstFaceIterator::ConstFaceIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Face::ConstFaceIterator::operator!=(const Face::ConstFaceIterator &it) const {
    return iter_ != it.iter_;
}

const Face &Face::ConstFaceIterator::operator*() const {
    return *iter_->rev()->face();
}

const Face *Face::ConstFaceIterator::ptr() const {
    return iter_->rev()->face();
}

const Face *Face::ConstFaceIterator::operator->() const {
    return iter_->rev()->face();
}

Face::ConstFaceIterator &Face::ConstFaceIterator::operator++() {
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Face::ConstFaceIterator Face::ConstFaceIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Face::ConstFaceIterator(tmp);
}

}  // namespace tinymesh
