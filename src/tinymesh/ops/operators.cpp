#define TINYMESH_API_EXPORT
#include "ops.h"

#include <vector>

#include "core/debug.h"
#include "core/vertex.h"
#include "core/halfedge.h"

namespace tinymesh {

void getMeshLaplacianAdjacent(const Mesh &mesh, EigenSparseMatrix &L) {
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        int valence = 0;
        for (auto it = v->v_begin(); it != v->v_end(); ++it) {
            const int j = it->index();
            triplets.emplace_back(i, j, -1.0);
            valence++;
        }
        triplets.emplace_back(i, i, valence);
    }
    L.resize(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());
}

void getMeshLaplacianCotangent(const Mesh &mesh, EigenSparseMatrix &L) {
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        double sumWgt = 0.0;
        for (auto it = v->he_begin(); it != v->he_end(); ++it) {
            const int j = it->dst()->index();
            const double weight = it->cotWeight();
            triplets.emplace_back(i, j, -weight);
            sumWgt += weight;
        }
        triplets.emplace_back(i, i, sumWgt);
    }
    L.resize(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());
}

/**
 * Reference:
 * Belkin et al., "Discrete Laplacian Operator on Meshed Surfaces," 2008.
 */
void getMeshLaplacianBelkin08(const Mesh &mesh, EigenSparseMatrix &L) {
    // Compute average edge length
    double avgEdgeLength = mesh.getMeanEdgeLength();

    // Construct sparse matrix
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    std::vector<EigenTriplet> areas;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        double sumWgt = 0.0;
        double area = 0.0;
        for (auto it = v->he_begin(); it != v->he_end(); ++it) {
            // Gaussian weight
            Vertex *u = it->dst();
            const double h = avgEdgeLength;
            const double norm = length(v->pos() - u->pos());
            const double weight = std::exp(-norm * norm / (4 * h)) / (4.0 * Pi * h);
            const int j = u->index();
            triplets.emplace_back(i, j, -weight);
            sumWgt += weight;
            // Area
            Vertex *w = it->next()->dst();
            area += length(cross(u->pos() - v->pos(), w->pos() - v->pos())) / 6.0;
        }
        triplets.emplace_back(i, i, sumWgt);
        areas.emplace_back(i, i, area);
    }
    EigenSparseMatrix W(N, N);
    W.setFromTriplets(triplets.begin(), triplets.end());
    L = W;

    // EigenSparseMatrix A(N, N);
    // A.setFromTriplets(areas.begin(), areas.end());
    // L = A * W;
}

EigenSparseMatrix getMeshLaplacian(const Mesh &mesh, MeshLaplace type) {
    EigenSparseMatrix L;
    if (type == MeshLaplace::Adjacent) {
        getMeshLaplacianAdjacent(mesh, L);
    } else if (type == MeshLaplace::Cotangent) {
        getMeshLaplacianCotangent(mesh, L);
    } else if (type == MeshLaplace::Belkin08) {
        getMeshLaplacianBelkin08(mesh, L);
    } else {
        Error("Unknown mesh Laplacian type specified: %d", (int)type);
    }
    return L;
}

}  // namespace tinymesh
