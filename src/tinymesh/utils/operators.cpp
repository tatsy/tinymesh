#define TINYMESH_API_EXPORT
#include "utils.h"

#include <vector>

#include "core/vertex.h"

namespace tinymesh {

void getMeshLaplacianAdjacent(Mesh &mesh, EigenSparseMatrix &L) {
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        for (auto it = v->v_begin(); it != v->v_end(); ++it) {
            const int j = it->index();
            triplets.emplace_back(i, j, 1.0);
        }
    }
    EigenSparseMatrix W(N, N);
    W.setFromTriplets(triplets.begin(), triplets.end());

    EigenSparseVector diags(N);
    for (int i = 0; i < N; i++) {
        diags += W.col(i);
    }

    W.diagonal() -= diags;
    // EigenSparseMatrix L = W - diags.asDiagonal();
}

void getMeshLaplacian(Mesh &mesh, MeshLaplace type, EigenSparseMatrix &L) {
    if (type == MeshLaplace::Adjacent) {
        getMeshLaplacianAdjacent(mesh, L);
    }
}

}  // namespace tinymesh
