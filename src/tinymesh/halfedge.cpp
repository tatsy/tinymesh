#define TINYMESH_API_EXPORT
#include "halfedge.h"

#include "vertex.h"

namespace tinymesh {

Halfedge::Halfedge() {}

Halfedge::Halfedge(const Halfedge &he) {
    this->operator=(he);
}

Halfedge::Halfedge(tinymesh::Halfedge &&he) noexcept {
    this->operator=(std::move(he));
}

double Halfedge::length() const {
    return (src()->pt() - dst()->pt()).length();
}


}  // namespace tinymesh
