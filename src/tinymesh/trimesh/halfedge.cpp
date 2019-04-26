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
    return ::length(src()->pt() - dst()->pt());
}


}  // namespace tinymesh
