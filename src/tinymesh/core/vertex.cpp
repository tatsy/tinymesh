#define TINYMESH_API_EXPORT
#include "vertex.h"

#include <set>
#include <vector>

#include "core/face.h"
#include "core/halfedge.h"

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

int Vertex::degree() const {
    int deg = 0;
    for (auto it = he_begin(); it != he_end(); ++it) {
        deg++;
    }
    return deg;
}

Vec3 Vertex::normal() const {
    Vec3 norm = Vec3(0.0);
    for (auto it = f_begin(); it != f_end(); ++it) {
        norm += it->normal() * it->area();
    }
    return normalize(norm);
}

bool Vertex::isBoundary() const {
    for (auto it = he_begin(); it != he_end(); ++it) {
        if (it->isBoundary()) {
            return true;
        }
    }
    return false;
}

double Vertex::K() const {
    std::vector<const Vertex *> neighbors;
    for (auto it = v_begin(); it != v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }

    const int N = static_cast<int>(neighbors.size());
    double sumAngles = 0.0;
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        const Vec3 e1 = neighbors[i]->pos() - pos();
        const Vec3 e2 = neighbors[j]->pos() - pos();
        sumAngles += std::atan2(length(cross(e1, e2)), dot(e1, e2));
        sumAreas += length(cross(e1, e2)) / 6.0;
    }
    return (2.0 * Pi - sumAngles) / sumAreas;
}

double Vertex::H() const {
    std::vector<const Vertex *> neighbors;
    for (auto it = v_begin(); it != v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }

    const int N = static_cast<int>(neighbors.size());
    Vec3 laplace = Vec3(0.0);
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int prev = (i - 1 + N) % N;
        const int post = (i + 1) % N;

        const Vec3 &p0 = pos();
        const Vec3 &p1 = neighbors[i]->pos();
        const Vec3 &p2 = neighbors[post]->pos();
        const Vec3 &p3 = neighbors[prev]->pos();

        const double sin_a = length(cross(p0 - p2, p1 - p2));
        const double cos_a = dot(p0 - p2, p1 - p2);
        const double cot_a = cos_a / std::max(sin_a, 1.0e-6);
        const double sin_b = length(cross(p0 - p3, p1 - p3));
        const double cos_b = dot(p0 - p3, p1 - p3);
        const double cot_b = cos_b / std::max(sin_b, 1.0e-6);
        const double weight = 0.5 * (cot_a + cot_b);
        laplace += weight * (neighbors[i]->pos() - pos());

        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        sumAreas += length(cross(e1, e2)) / 6.0;
    }
    return -1.0 * dot(laplace, this->normal()) / (2.0 * sumAreas);
}

Vertex::VertexIterator Vertex::v_begin() {
    return Vertex::VertexIterator(halfedge_);
}

Vertex::VertexIterator Vertex::v_end() {
    return Vertex::VertexIterator(nullptr);
}

Vertex::ConstVertexIterator Vertex::v_begin() const {
    return Vertex::ConstVertexIterator(halfedge_);
}

Vertex::ConstVertexIterator Vertex::v_end() const {
    return Vertex::ConstVertexIterator(nullptr);
}

Vertex::HalfedgeIterator Vertex::he_begin() {
    return Vertex::HalfedgeIterator(halfedge_);
}

Vertex::HalfedgeIterator Vertex::he_end() {
    return Vertex::HalfedgeIterator(nullptr);
}

Vertex::ConstHalfedgeIterator Vertex::he_begin() const {
    return Vertex::ConstHalfedgeIterator(halfedge_);
}

Vertex::ConstHalfedgeIterator Vertex::he_end() const {
    return Vertex::ConstHalfedgeIterator(nullptr);
}

Vertex::FaceIterator Vertex::f_begin() {
    return Vertex::FaceIterator(halfedge_);
}

Vertex::FaceIterator Vertex::f_end() {
    return Vertex::FaceIterator(nullptr);
}

Vertex::ConstFaceIterator Vertex::f_begin() const {
    return Vertex::ConstFaceIterator(halfedge_);
}

Vertex::ConstFaceIterator Vertex::f_end() const {
    return Vertex::ConstFaceIterator(nullptr);
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

Vertex *Vertex::VertexIterator::ptr() const {
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
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::VertexIterator(tmp);
}

// ----------
// ConstVertexIterator
// ----------

Vertex::ConstVertexIterator::ConstVertexIterator(tinymesh::Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::ConstVertexIterator::operator!=(const Vertex::ConstVertexIterator &it) const {
    return iter_ != it.iter_;
}

const Vertex &Vertex::ConstVertexIterator::operator*() const {
    return *iter_->dst();
}

const Vertex *Vertex::ConstVertexIterator::ptr() const {
    return iter_->dst();
}

const Vertex *Vertex::ConstVertexIterator::operator->() const {
    return iter_->dst();
}

Vertex::ConstVertexIterator &Vertex::ConstVertexIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::ConstVertexIterator Vertex::ConstVertexIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::ConstVertexIterator(tmp);
}

// ----------
// HalfedgeIterator
// ----------

Vertex::HalfedgeIterator::HalfedgeIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::HalfedgeIterator::operator!=(const Vertex::HalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

Halfedge &Vertex::HalfedgeIterator::operator*() {
    return *iter_;
}

Halfedge *Vertex::HalfedgeIterator::ptr() const {
    return iter_;
}

Halfedge *Vertex::HalfedgeIterator::operator->() const {
    return iter_;
}

Vertex::HalfedgeIterator &Vertex::HalfedgeIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::HalfedgeIterator Vertex::HalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::HalfedgeIterator(tmp);
}

// ----------
// ConstHalfedgeIterator
// ----------

Vertex::ConstHalfedgeIterator::ConstHalfedgeIterator(Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::ConstHalfedgeIterator::operator!=(const Vertex::ConstHalfedgeIterator &it) const {
    return iter_ != it.iter_;
}

const Halfedge &Vertex::ConstHalfedgeIterator::operator*() const {
    return *iter_;
}

const Halfedge *Vertex::ConstHalfedgeIterator::ptr() const {
    return iter_;
}

const Halfedge *Vertex::ConstHalfedgeIterator::operator->() const {
    return iter_;
}

Vertex::ConstHalfedgeIterator &Vertex::ConstHalfedgeIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::ConstHalfedgeIterator Vertex::ConstHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::ConstHalfedgeIterator(tmp);
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

// ----------
// ConstFaceIterator
// ----------

Vertex::ConstFaceIterator::ConstFaceIterator(tinymesh::Halfedge *he)
    : halfedge_{ he }
    , iter_{ he } {
}

bool Vertex::ConstFaceIterator::operator!=(const Vertex::ConstFaceIterator &it) const {
    return iter_ != it.iter_;
}

const Face &Vertex::ConstFaceIterator::operator*() const {
    return *iter_->face();
}

const Face *Vertex::ConstFaceIterator::ptr() const {
    return iter_->face();
}

const Face *Vertex::ConstFaceIterator::operator->() const {
    return iter_->face();
}

Vertex::ConstFaceIterator &Vertex::ConstFaceIterator::operator++() {
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::ConstFaceIterator Vertex::ConstFaceIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->next();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::ConstFaceIterator(tmp);
}

}  // namespace tinymesh
