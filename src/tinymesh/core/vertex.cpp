#define TINYMESH_API_EXPORT
#include "vertex.h"

#include <set>
#include <vector>

#include "core/debug.h"
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
    // Normal weighting scheme in [Max 1999]
    // Reference: Weights for Computing Vertex Normals from Facet Normals
    Vec3 norm = Vec3(0.0);
    for (auto it = he_begin(); it != he_end(); ++it) {
        const Vec3 p0 = this->pos();
        const Vec3 p1 = it->dst()->pos();
        const Vec3 p2 = it->next()->dst()->pos();
        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        const double l1 = length(e1);
        const double l2 = length(e2);
        norm += cross(e1, e2) / (l1 * l1 * l2 * l2);
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

double Vertex::volonoiArea() const {
    std::vector<const Vertex *> neighbors;
    for (auto it = v_begin(); it != v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }

    const int N = (int)neighbors.size();
    double sumArea = 0.0;
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        const Vec3 &p0 = pos();
        const Vec3 &p1 = neighbors[i]->pos();
        const Vec3 &p2 = neighbors[j]->pos();
        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;

        // Safe Volonoi area calculation [Meyer et al. 2003]
        if (!obtuse(p0, p1, p2)) {
            const double l1 = length(e1);
            const double l2 = length(e2);
            const double A1 = 0.125 * l1 * l1 * cot(p0, p2, p1);
            const double A2 = 0.125 * l2 * l2 * cot(p0, p1, p2);
            sumArea += A1 + A2;
        } else {
            if (dot(e1, e2) < 0.0) {
                // triangle is obtuse at "p0"
                sumArea += 0.5 * length(cross(e1, e2)) * 0.5;
            } else {
                // otherwise
                sumArea += 0.5 * length(cross(e1, e2)) * 0.25;
            }
        }
    }
    return sumArea;
}

double Vertex::volonoiArea(const Face *const f) const {
    int index = -1;
    int count = 0;

    std::vector<const Vertex *> vertices;
    for (auto it = f->v_begin(); it != f->v_end(); ++it) {
        vertices.push_back(it.ptr());
        if (it.ptr() == this) {
            index = count;
        }
        count += 1;
    }
    Assertion(index >= 0, "Vertex is not included in the face!");
    Assertion(vertices.size() == 3, "Non-triangle face detected!");

    const Vec3 p0 = vertices[index]->pos();
    const Vec3 p1 = vertices[(index + 1) % vertices.size()]->pos();
    const Vec3 p2 = vertices[(index + 2) % vertices.size()]->pos();
    const Vec3 e1 = p1 - p0;
    const Vec3 e2 = p2 - p0;

    // Safe Volonoi area calculation [Meyer et al. 2003]
    double area = 0.0;
    if (!obtuse(p0, p1, p2)) {
        const double l1 = length(e1);
        const double l2 = length(e2);
        const double A1 = 0.125 * l1 * l1 * cot(p0, p2, p1);
        const double A2 = 0.125 * l2 * l2 * cot(p0, p1, p2);
        area = A1 + A2;
    } else {
        if (dot(e1, e2) < 0.0) {
            // triangle is obtuse at "p0"
            area = 0.5 * length(cross(e1, e2)) * 0.5;
        } else {
            // otherwise
            area = 0.5 * length(cross(e1, e2)) * 0.25;
        }
    }
    return area;
}

double Vertex::K() const {
    // NOTE: VertexIterator traverses neighboring vertices in the clockwise order.
    std::vector<const Vertex *> neighbors;
    for (auto it = v_begin(); it != v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }
    std::reverse(neighbors.begin(), neighbors.end());

    const int N = static_cast<int>(neighbors.size());
    double sumAngles = 0.0;
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int j = (i + 1) % N;
        const Vec3 &p0 = pos();
        const Vec3 &p1 = neighbors[i]->pos();
        const Vec3 &p2 = neighbors[j]->pos();
        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        sumAngles += std::atan2(length(cross(e1, e2)), dot(e1, e2));

        // Safe Volonoi area calculation [Meyer et al. 2003]
        if (!obtuse(p0, p1, p2)) {
            const double l1 = length(e1);
            const double l2 = length(e2);
            const double A1 = 0.125 * l1 * l1 * cot(p0, p2, p1);
            const double A2 = 0.125 * l2 * l2 * cot(p0, p1, p2);
            sumAreas += A1 + A2;
        } else {
            if (dot(e1, e2) < 0.0) {
                // triangle is obtuse at "p0"
                sumAreas += 0.5 * length(cross(e1, e2)) * 0.5;
            } else {
                // otherwise
                sumAreas += 0.5 * length(cross(e1, e2)) * 0.25;
            }
        }
    }
    return (2.0 * Pi - sumAngles) / sumAreas;
}

double Vertex::H() const {
    // NOTE: VertexIterator traverses neighboring vertices in the clockwise order.
    std::vector<const Vertex *> neighbors;
    for (auto it = v_begin(); it != v_end(); ++it) {
        neighbors.push_back(it.ptr());
    }
    std::reverse(neighbors.begin(), neighbors.end());

    const int N = static_cast<int>(neighbors.size());
    Vec3 laplace = Vec3(0.0);
    double sumAreas = 0.0;
    for (int i = 0; i < N; i++) {
        const int prev = (i - 1 + N) % N;
        const int next = (i + 1) % N;

        const Vec3 &p0 = pos();
        const Vec3 &p1 = neighbors[i]->pos();
        const Vec3 &p2 = neighbors[next]->pos();
        const Vec3 &p3 = neighbors[prev]->pos();

        const double cot_a = cot(p0, p2, p1);
        const double cot_b = cot(p0, p3, p1);
        const double weight = cot_a + cot_b;
        laplace += weight * (neighbors[i]->pos() - pos());

        // Safe Volonoi area calculation [Meyer et al. 2003]
        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        if (!obtuse(p0, p1, p2)) {
            const double l1 = length(e1);
            const double l2 = length(e2);
            const double A1 = 0.125 * l1 * l1 * cot(p0, p2, p1);
            const double A2 = 0.125 * l2 * l2 * cot(p0, p1, p2);
            sumAreas += A1 + A2;
        } else {
            if (dot(e1, e2) < 0.0) {
                // triangle is obtuse at "p0"
                sumAreas += 0.5 * length(cross(e1, e2)) * 0.5;
            } else {
                // otherwise
                sumAreas += 0.5 * length(cross(e1, e2)) * 0.25;
            }
        }
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
