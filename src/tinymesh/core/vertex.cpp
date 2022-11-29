#define TINYMESH_API_EXPORT
#include "vertex.h"

#include <set>
#include <vector>

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
        deg++;
    }
    return deg;
}

Vec3 Vertex::normal() {
    Vec3 norm = Vec3(0.0);
    for (auto it = f_begin(); it != f_end(); ++it) {
        norm += it->normal() * it->area();
    }
    return normalize(norm);
}

double Vertex::K() {
    std::vector<Vertex *> neighbors;
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

double Vertex::H() {
    std::vector<Vertex *> neighbors;
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

Vertex *Vertex::VertexIterator::ptr() const {
    return iter_->dst();
}

Vertex *Vertex::VertexIterator::operator->() const {
    return iter_->dst();
}

Vertex::VertexIterator &Vertex::VertexIterator::operator++() {
    iter_ = iter_->prev()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::VertexIterator Vertex::VertexIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->prev()->rev();
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
    iter_ = iter_->rev()->prev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::InHalfedgeIterator Vertex::InHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->rev()->prev();
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
    iter_ = iter_->prev()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::OutHalfedgeIterator Vertex::OutHalfedgeIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->prev()->rev();
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
    iter_ = iter_->prev()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return *this;
}

Vertex::FaceIterator Vertex::FaceIterator::operator++(int) {
    Halfedge *tmp = iter_;
    iter_ = iter_->prev()->rev();
    if (iter_ == halfedge_) {
        iter_ = nullptr;
    }
    return Vertex::FaceIterator(tmp);
}

}  // namespace tinymesh
