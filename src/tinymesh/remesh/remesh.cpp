#define TINYMESH_API_EXPORT
#include "remesh.h"

#include <random>
#include <algorithm>
#include <atomic>

#include "core/vec.h"
#include "core/debug.h"
#include "core/openmp.h"
#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"
#include "filters/filters.h"

namespace tinymesh {

void remeshTriangular(Mesh &mesh, double shortLength, double longLength, double keepAngleLessThan, int iterations, bool verbose) {
    Assertion(mesh.verify(), "Invalid mesh!");

    // Compute average edge length
    int count;
    double Lavg, Lvar;
    Lavg = 0.0;
    Lvar = 0.0;

    count = 0;
    for (int i = 0; i < (int)mesh.numHalfedges(); i++) {
        Halfedge *he = mesh.halfedge(i);
        const double l = he->length();
        Lavg += l;
        Lvar += l * l;
        count += 1;
    }

    Lavg = Lavg / count;
    Lvar = Lvar / count - Lavg * Lavg;

    // Check whether each vertex is on feature line
    std::vector<double> minDiheds(mesh.numVertices(), Pi);
    for (int i = 0; i < (int)mesh.numVertices(); i++) {
        Vertex *v = mesh.vertex(i);
        std::vector<Vec3> neighbors;
        for (auto vit = v->v_begin(); vit != v->v_end(); ++vit) {
            neighbors.push_back(vit->pos());
        }

        const auto nn = static_cast<int>(neighbors.size());
        for (int j = 0; j < nn; j++) {
            const int k = (j + 1) % nn;
            const int l = (j - 1 + nn) % nn;

            const Vec3 p0 = v->pos();
            const Vec3 p1 = neighbors[j];
            const Vec3 p2 = neighbors[k];
            const Vec3 p3 = neighbors[l];
            const double dihed = dihedral(p2, p0, p1, p3);
            minDiheds[i] = std::min(dihed, minDiheds[i]);
        }

        if (minDiheds[i] < keepAngleLessThan) {
            v->setIsStatic(true);
        }
    }

    // Initialize random number generator
    std::vector<uint32_t> indices;
    std::random_device randev;
    std::mt19937 rnd(randev());

    // Remesh loop
    for (int k = 0; k < iterations; k++) {
        if (verbose) {
            Info("*** Original #%d ***\n", k + 1);
            Info("#vert: %d\n", (int)mesh.numVertices());
            Info("#face: %d\n", (int)mesh.numFaces());
        }

        // Split long edges
        indices.clear();
        for (int i = 0; i < (int)mesh.numHalfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < (int)mesh.numHalfedges()) {
                Halfedge *he = mesh.halfedge(i);
                const Vec3 p1 = he->src()->pos();
                const Vec3 p2 = he->dst()->pos();
                const double l = length(p1 - p2);

                if (l >= Lavg * longLength) {
                    mesh.splitHE(he);
                }
            }
        }

        if (verbose) {
            Info("*** After split ***\n");
            Info("#vert: %d\n", (int)mesh.numVertices());
            Info("#face: %d\n", (int)mesh.numFaces());
        }

        mesh.verify();

        // Collapse short edges
        indices.clear();
        for (int i = 0; i < (int)mesh.numHalfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < (int)mesh.numHalfedges()) {
                Halfedge *he = mesh.halfedge(i);
                if (he->face()->isStatic() || he->rev()->face()->isStatic()) {
                    continue;
                }

                const int i1 = he->src()->index();
                const int i2 = he->dst()->index();
                if (he->src()->isStatic() || he->dst()->isStatic()) {
                    continue;
                }

                const Vec3 p1 = he->src()->pos();
                const Vec3 p2 = he->dst()->pos();
                const double l = length(p1 - p2);

                if (l <= Lavg * shortLength) {
                    // Check if collapse does not generate long edge
                    Vertex *a = he->src();
                    Vertex *b = he->dst();
                    bool collapseOK = true;
                    for (auto vit = b->v_begin(); vit != b->v_end(); ++vit) {
                        if (length(a->pos() - vit->pos()) >= Lavg * longLength) {
                            collapseOK = false;
                            break;
                        }
                    }

                    // Collapse
                    if (collapseOK) {
                        mesh.collapseHE(he);
                    }
                }
            }
        }

        if (verbose) {
            Info("*** After collapse ***\n");
            Info("#vert: %d\n", (int)mesh.numVertices());
            Info("#face: %d\n", (int)mesh.numFaces());
        }

        // Flip edges
        for (int i = 0; i < (int)mesh.numHalfedges(); i++) {
            Halfedge *he = mesh.halfedge(i);
            if (he->face()->isBoundary() || he->rev()->face()->isBoundary()) {
                continue;
            }

            if (he->face()->isStatic() || he->rev()->face()->isStatic()) {
                continue;
            }

            Vertex *v0 = he->src();
            Vertex *v1 = he->dst();
            Vertex *v2 = he->next()->dst();
            Vertex *v3 = he->rev()->next()->dst();

            if (v0->isStatic() || v1->isStatic()) {
                continue;
            }

            const int d0 = v0->degree();
            const int d1 = v1->degree();
            const int d2 = v2->degree();
            const int d3 = v3->degree();
            const int t0 = v0->isBoundary() ? 4 : 6;
            const int t1 = v1->isBoundary() ? 4 : 6;
            const int t2 = v2->isBoundary() ? 4 : 6;
            const int t3 = v3->isBoundary() ? 4 : 6;

            const int score = std::abs(d0 - t0) + std::abs(d1 - t1) + std::abs(d2 - t2) + std::abs(d3 - t3);
            const int after = std::abs(d0 - 1 - t0) + std::abs(d1 - 1 - t1) + std::abs(d2 + 1 - t2) + std::abs(d3 + 1 - t3);
            if (score > after) {
                mesh.flipHE(he);
            }
        }

        // Smoothing
        smoothTaubin(mesh);
    }
}

}  // namespace tinymesh