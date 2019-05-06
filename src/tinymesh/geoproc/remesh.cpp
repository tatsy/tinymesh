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

void remeshIncremental(Mesh &mesh, double ratioLower, double ratioUpper, int maxiter) {
    Assertion(mesh.verify(), "Invalid mesh!");

    // Compute average edge length
    int count;
    double Lavg, Lvar, Lstd;
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
    Lstd = std::sqrt(Lvar);

    // Initialize random number generator
    std::vector<uint32_t> indices;
    std::random_device randev;
    std::mt19937 rnd(randev());

    // Remesh loop
    for (int k = 0; k < maxiter; k++) {
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
                const double l = he->length();
                const double p = (l - Lavg) / Lstd;
                if (l >= Lavg * ratioUpper) {
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
                const double l = he->length();
                const double p = (l - Lavg) / Lstd;
                if (l <= Lavg * ratioLower) {
                    // Check if collapse does not generate long edge
                    Vertex *a = he->src();
                    Vertex *b = he->dst();
                    bool collapseOK = true;
                    for (auto vit = b->v_begin(); vit != b->v_end(); ++vit) {
                        if (length(a->pos() - vit->pos()) >= Lavg * ratioUpper) {
                            collapseOK = false;
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

            Vertex *v0 = he->src();
            Vertex *v1 = he->dst();
            Vertex *v2 = he->next()->dst();
            Vertex *v3 = he->rev()->next()->dst();
            const int d0 = v0->degree();
            const int d1 = v1->degree();
            const int d2 = v2->degree();
            const int d3 = v3->degree();

            const int score = std::abs(d0 - 6) + std::abs(d1 - 6) + std::abs(d2 - 6) + std::abs(d3 - 6);
            const int after = std::abs(d0 - 1 - 6) + std::abs(d1 - 1 - 6) + std::abs(d2 + 1 - 6) + std::abs(d3 + 1 - 6);
            if (score > after) {
                mesh.flipHE(he);
            }
        }

        // Smoothing
        for (int loop = 0; loop < 3; loop++) {
            smooth(mesh);
        }
    }
}

}  // namespace tinymesh