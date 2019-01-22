#define TINYMESH_API_EXPORT
#include "halfedge.h"

namespace tinymesh {

Halfedge::Halfedge() {}

Halfedge::Halfedge(const Halfedge &he) {
    this->operator=(he);
}

Halfedge::Halfedge(tinymesh::Halfedge &&he) noexcept {
    this->operator=(std::move(he));
}

}  // namespace tinymesh
