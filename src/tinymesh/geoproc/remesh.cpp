#define TINYMESH_API_EXPORT
#include "remesh.h"

#include <random>
#include <algorithm>
#include <atomic>

#include "core/vec.h"
#include "core/openmp.h"
#include "polymesh/mesh.h"
#include "polymesh/vertex.h"
#include "polymesh/halfedge.h"
#include "polymesh/face.h"
#include "geoproc/smooth.h"

namespace tinymesh {

void remeshIncremental(Mesh &mesh, double shortLength, double longLength, double angleThresh, int iterations) {
    Assertion(mesh.verify(), "Invalid mesh!");

    // Compute average edge length
    int count;
    double Lavg, Lvar;
    Lavg = 0.0;
    Lvar = 0.0;

    count = 0;
    for (int i = 0; i < mesh.num_halfedges(); i++) {
        Halfedge *he = mesh.halfedge(i);
        const double l = he->length();
        Lavg += l;
        Lvar += l * l;
        count += 1;
    }

    Lavg = Lavg / count;
    Lvar = Lvar / count - Lavg * Lavg;

    // Initialize random number generator
    std::vector<uint32_t> indices;
    std::random_device randev;
    std::mt19937 rnd(randev());

    // Remesh loop
    for (int k = 0; k < iterations; k++) {
        printf("*** Original #%d ***\n", k + 1);
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Split long edges
        indices.clear();
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < mesh.num_halfedges()) {
                Halfedge *he = mesh.halfedge(i);
                if (he->face()->isStatic() || he->rev()->face()->isStatic()) {
                    continue;
                }

                const Vec3 p1 = he->src()->pos();
                const Vec3 p2 = he->dst()->pos();
                const double l = length(p1 - p2);

                if (l >= Lavg * longLength) {
                    mesh.splitHE(he);
                }
            }
        }

        printf("*** After split ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        mesh.verify();

        // Collapse short edges
        indices.clear();
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < mesh.num_halfedges()) {
                Halfedge *he = mesh.halfedge(i);
                if (he->face()->isStatic() || he->rev()->face()->isStatic()) {
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

        printf("*** After collapse ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Flip edges
        for (int i = 0; i < mesh.num_halfedges(); i++) {
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

            // Compute dihedral angle before and after flip
            const Vec3 p0 = v0->pos();
            const Vec3 p1 = v1->pos();
            const Vec3 p2 = v2->pos();
            const Vec3 p3 = v3->pos();
            const Vec3 n0 = cross(p1 - p0, p2 - p0);
            const Vec3 n1 = cross(p1 - p0, p3 - p0);
            if (length(n0) == 0.0 || length(n1) == 0.0) {
                continue;
            }

            if (dot(n0, n1) / (length(n0) * length(n1)) > angleThresh) continue;

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
        laplace_smooth(mesh);
    }
}

}  // namespace tinymesh