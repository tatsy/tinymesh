#define TINYMESH_API_EXPORT
#include "face.h"

#include "halfedge.h"

namespace tinymesh {

Face::Face() {   
}


Face::VertexIterator Face::v_begin() {
    return Face::VertexIterator(m_he);
}

Face::VertexIterator Face::v_end() {
    return Face::VertexIterator(nullptr);
}

// ----------
// VertexIterator
// ----------

Face::VertexIterator::VertexIterator(Halfedge *he)
    : m_he{ he }
    , m_init{ he } {
}

bool Face::VertexIterator::operator!=(const Face::VertexIterator &it) const {
    return m_he != it.m_he;
}

Vertex &Face::VertexIterator::operator*() {
    return *m_he->dst();
}

Vertex *Face::VertexIterator::ptr() const {
    return m_he->dst();
}

Vertex *Face::VertexIterator::operator->() const {
    return m_he->dst();
}

Face::VertexIterator &Face::VertexIterator::operator++() {
    m_he = m_he->rev()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return *this;
}

Face::VertexIterator Face::VertexIterator::operator++(int) {
    Halfedge *tmp = m_he;
    m_he = m_he->rev()->next();
    if (m_he == m_init) {
        m_he = nullptr;
    }
    return Face::VertexIterator(tmp);
}

}  // namespace tinymesh
