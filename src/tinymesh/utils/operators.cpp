#define TINYMESH_API_EXPORT
#include "utils.h"

#include <vector>

#include "core/debug.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/halfedge.h"

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
    L.resize(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());

    EigenSparseVector diags(N);
    for (int i = 0; i < N; i++) {
        diags += L.col(i);
    }

    L.diagonal() -= diags;
}

void getMeshLaplacianCotangent(Mesh &mesh, EigenSparseMatrix &L) {
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        for (auto it = v->ohe_begin(); it != v->ohe_end(); ++it) {
            Halfedge &he = *it;
            Vertex *u = he.dst();
            Vertex *l = he.next()->dst();
            Vertex *r = he.rev()->next()->dst();

            const Vec3 &p0 = v->pos();
            const Vec3 &p1 = u->pos();
            const Vec3 &p2 = l->pos();
            const Vec3 &p3 = r->pos();
            const double sin_a = length(cross(p0 - p2, p1 - p2));
            const double cos_a = dot(p0 - p2, p1 - p2);
            const double cot_a = cos_a / std::max(sin_a, 1.0e-6);
            const double sin_b = length(cross(p0 - p3, p1 - p3));
            const double cos_b = dot(p0 - p3, p1 - p3);
            const double cot_b = cos_b / std::max(sin_b, 1.0e-6);
            const double weight = 0.5 * (cot_a + cot_b);

            const int j = u->index();
            triplets.emplace_back(i, j, weight);
        }
    }
    L.resize(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());

    EigenSparseVector diags(N);
    for (int i = 0; i < N; i++) {
        diags += L.col(i);
    }

    L.diagonal() -= diags;
}

/**
 * Reference:
 * Belkin et al., "Discrete Laplacian Operator on Meshed Surfaces," 2008.
 */
void getMeshLaplacianBelkin08(Mesh &mesh, EigenSparseMatrix &L) {
    // Compute average edge length
    double avgEdgeLength = 0.0;
    for (int i = 0; i < mesh.numEdges(); i++) {
        avgEdgeLength += mesh.edge(i)->length();
    }
    avgEdgeLength /= mesh.numEdges();

    // Construct sparse matrix
    const int N = mesh.numVertices();
    std::vector<EigenTriplet> triplets;
    std::vector<EigenTriplet> areas;
    for (int i = 0; i < N; i++) {
        Vertex *v = mesh.vertex(i);
        double area = 0.0;
        for (auto it = v->ohe_begin(); it != v->ohe_end(); ++it) {
            // Gaussian weight
            Vertex *u = it->dst();
            const double h = avgEdgeLength;
            const double norm = length(v->pos() - u->pos());
            const double weight = std::exp(-norm * norm / (4 * h)) / std::pow(4.0 * Pi * h, 1.5);
            const int j = it->index();
            triplets.emplace_back(i, j, weight);
            // Area
            Vertex *w = it->next()->dst();
            area += 0.5 * length(cross(u->pos() - v->pos(), w->pos() - v->pos())) / 3.0;
        }
        areas.emplace_back(i, i, area);
    }
    L.resize(N, N);
    L.setFromTriplets(triplets.begin(), triplets.end());

    EigenSparseVector diags(N);
    for (int i = 0; i < N; i++) {
        diags += L.col(i);
    }
    L.diagonal() -= diags;

    EigenSparseMatrix A(N, N);
    A.setFromTriplets(areas.begin(), areas.end());
    L *= A;
}

EigenSparseMatrix getMeshLaplacian(Mesh &mesh, MeshLaplace type) {
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
