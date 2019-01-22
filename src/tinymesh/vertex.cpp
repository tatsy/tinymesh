#define TINYMESH_API_EXPORT
#include "vertex.h"

#include "halfedge.h"

namespace tinymesh {

// ----------
// Vertex
// ----------

Vertex::Vertex() {}
Vertex::Vertex(const Vector &v) : m_pt{ v } {}

Vertex::VertexIterator Vertex::v_begin() {
    return Vertex::VertexIterator(m_he);
}

Vertex::VertexIterator Vertex::v_end() {
    return Vertex::VertexIterator(nullptr);
}

Vertex::InHalfedgeIterator Vertex::ihe_begin() {
    return Vertex::InHalfedgeIterator(m_he);
}

Vertex::InHalfedgeIterator Vertex::ihe_end() {
    return Vertex::InHalfedgeIterator(nullptr);
}

Vertex::OutHalfedgeIterator Vertex::ohe_begin() {
    return Vertex::OutHalfedgeIterator(m_he);
}

Vertex::OutHalfedgeIterator Vertex::ohe_end() {
    return Vertex::OutHalfedgeIterator(nullptr);
}

// ----------
// VertexIterator
// ----------

Vertex::VertexIterator::VertexIterator(tinymesh::Halfedge *he)
    : m_he{ he }
    , m_init{ he } {
}

bool Vertex::VertexIterator::operator!=(const Vertex::VertexIterator &it) const {
    return m_he != it.m_he;
}

Vertex &Vertex::VertexIterator::operator*() {
    return *m_he->v_end();
}

Vertex *Vertex::VertexIterator::operator->() const {
    return m_he->v_end();
}

Vertex::VertexIterator &Vertex::VertexIterator::operator++() {
    m_he = m_he->opposite()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return *this;
}

Vertex::VertexIterator Vertex::VertexIterator::operator++(int) {
    Halfedge *tmp = m_he;
    m_he = m_he->opposite()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return Vertex::VertexIterator(tmp);
}

// ----------
// InHalfedgeIterator
// ----------

Vertex::InHalfedgeIterator::InHalfedgeIterator(tinymesh::Halfedge *he)
    : m_he{ he }
    , m_init{ he } {
}

bool Vertex::InHalfedgeIterator::operator!=(const Vertex::InHalfedgeIterator &it) const {
    return m_he != it.m_he;
}

Halfedge &Vertex::InHalfedgeIterator::operator*() {
    return *m_he;
}

Halfedge *Vertex::InHalfedgeIterator::operator->() const {
    return m_he;
}

Vertex::InHalfedgeIterator &Vertex::InHalfedgeIterator::operator++() {
    m_he = m_he->next()->opposite();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return *this;
}

Vertex::InHalfedgeIterator Vertex::InHalfedgeIterator::operator++(int) {
    Halfedge *tmp = m_he;
    m_he = m_he->next()->opposite();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return Vertex::InHalfedgeIterator(tmp);
}

// ----------
// OutHalfedgeIterator
// ----------

Vertex::OutHalfedgeIterator::OutHalfedgeIterator(Halfedge *he)
    : m_he{ he }
    , m_init{ he } {
}

bool Vertex::OutHalfedgeIterator::operator!=(const Vertex::OutHalfedgeIterator &it) const {
    return m_he != it.m_he;
}

Halfedge &Vertex::OutHalfedgeIterator::operator*() {
    return *m_he;
}

Halfedge *Vertex::OutHalfedgeIterator::operator->() const {
    return m_he;
}

Vertex::OutHalfedgeIterator &Vertex::OutHalfedgeIterator::operator++() {
    m_he = m_he->opposite()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return *this;
}

Vertex::OutHalfedgeIterator Vertex::OutHalfedgeIterator::operator++(int) {
    Halfedge *tmp = m_he;
    m_he = m_he->opposite()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return Vertex::OutHalfedgeIterator(tmp);
}

}  // naemspace tinymesh
